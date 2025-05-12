################################# #
# maxnet ENMdetails object ####
################################# #

# This function predicts with a cross-validated Maxent model
predict.CVmaxnet <- function(
    object, newdata, clamp = TRUE,
    type = c("link", "exponential", "cloglog", "logistic"),
    offset = NULL, ...) {
  # |GEPB|: First lines in predict.maxnet()
  # na_action <- options("na.action")[[1]]
  # on.exit(options(na.action = na_action))
  # options(na.action = "na.pass")
  
  newdata <- data.frame(newdata) # sparse matrices cause issues
  if (clamp) {
    for (v in intersect(names(object$varmax), names(newdata))) {
      newdata[, v] <- pmin(pmax(newdata[, v], object$varmin[v]),
                           object$varmax[v])
    }
  }
  
  # maxnet has a bug where some hinges are duplicated, so combining them
  #  cleaning up the formatting of these types for use with formula
  b.names <- names(object$betas)
  to.combine <- b.names[which(duplicated(b.names))]
  # CHECK GEPB betas
  if (length(to.combine) > 0) {
    for (tc in to.combine) {
      tmp <- which(b.names %in% tc)
      object$betas[tmp] <- sum(object$betas[tmp])
      object$betas <- object$betas[-tmp[-1]]
      b.names <- names(object$betas)
    }
  }
  if (!length(unique(names(object$betas))) == length((names(object$betas)))) {
    print('some duplicated feature names still')
  }
  
  terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)",
               names(object$betas))
  terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\\2)",
               terms)
  terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)",
               terms)
  # changed f to form because f is a command when using browser()
  form <- formula(paste("~", paste(terms, collapse = " + "), "-1"))
  
  mm <- model.matrix(form, data.frame(newdata))
  if (clamp) {
    mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                 object$featuremaxs[names(object$betas)]))
  }
  # CM added offset; CM already logs this on input, so this isn't general
  # CM ditched alpha because it required a wasted projection step in maxnet, and
  # its just a normalizing constant that you can take care of later, i think.
  # maybe its an issue for the rescaled ones, but don't care...
  t1 <- proc.time()
  link <- (mm %*% object$betas) + ifelse(!is.null(offset), offset, 0) + object$alpha
  t2 <-  proc.time() - t1
  type <- match.arg(type)
  if (type == "link") {
    return(link)
  }
  if (type == "exponential") {
    out <- exp(link)
    return(out)
  }
  if (type == "cloglog") {
    return(1 - exp(0 - exp(object$entropy + link)))
  }
  if (type == "logistic") {
    return(1 / (1 + exp(-object$entropy - link)))
  }
}

categoricalval <- function(x, category) {
  ifelse(x == category, 1, 0)
}

thresholdval <- function(x, knot) {
  ifelse(x >= knot, 1, 0)
}

hingeval <- function(x, min, max) {
  pmin(1, pmax(0, (x - min) / (max - min)))
}

# This function fits a cross-validated Maxent model
CVmaxnet <- function(p, data, f = maxnet::maxnet.formula(p, data), regmult = 1,
                     regfun = maxnet::maxnet.default.regularization,
                     weights = p + (1 - p) * 100, offset = NULL, maxit = 2e5,
                     standardize = TRUE, nfolds = 10, foldid = NULL,
                     parallel = FALSE, verbose = TRUE, predictions = TRUE,
                     quadraticUpperLimit0 = FALSE, lambda.sel = "1se", ...) {
  
  # Remove NAs of data and p
  keep <- which(complete.cases(data))
  data <- data[keep, ]
  p <- p[keep]
  weights <- weights[keep]
  offset <- offset[keep]
  
  mm <- model.matrix(f, data)
  reg <- regfun(p = p, mm) * regmult
  glmnet::glmnet.control(pmin = 1e-08,
                         fdev = 1e-5)
  lambda <- 10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*sum(p)/sum(weights)
  
  # set upper.limits on quadratic coeffs to force them to be negative
  if (quadraticUpperLimit0) {
    cmquadraticUpperLimit0 <- function(m) {
      up <- rep(Inf, ncol(m))
      isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
      to0 <- isquadratic(colnames(m))
      up[to0] <- 0
      up
    }
    upper.limits <-  cmquadraticUpperLimit0(mm)
  } else {
    upper.limits <-  rep(Inf, ncol(mm))
  }
  
  if (verbose) print('Starting cross-validated model.')
  
  # Fit CV glmnet model
  model.cv <- glmnet::cv.glmnet(x = mm, y = p, family = "binomial",
                        standardize = standardize, penalty.factor = reg,
                        nlambda = 200, lambda = lambda,
                        weights = weights, offset = offset, maxit = maxit,
                        nfolds = nfolds, foldid = foldid, parallel = FALSE,
                        upper.limits = upper.limits, type.measure = 'deviance', ...)
  class(model.cv) <- c("CVmaxnet", class(model.cv))
  
  # SELECT LAMBDA
  # Rules for choosing lambda.
  lam.1se.index <- which(model.cv$lambda == model.cv$lambda.1se)
  lam.min.index <- which(model.cv$lambda == model.cv$lambda.min)
  ## Calculate  0.5 se and 0.25 se manually
  se <- model.cv$cvsd[lam.min.index]
  min.cvm <- model.cv$cvm[lam.min.index]
  # 0.5 se
  cvm.cut.half.se <- min.cvm + .5 * se
  which.cut.rev.half.se <- which.min(abs(model.cv$cvm - (min.cvm + se * 0.5)))
  lambda.val.half.se <- model.cv$lambda[which.cut.rev.half.se]
  lambda.0.5se.index <- which(model.cv$lambda == lambda.val.half.se)
  model.cv$index[3] <- lambda.0.5se.index
  # 0.25 se
  cvm.cut.quar.se <- min.cvm + .25 * se
  which.cut.rev.quar.se <- which.min(abs(model.cv$cvm - (min.cvm + se * 0.25)))
  lambda.val.quar.se <- model.cv$lambda[which.cut.rev.quar.se]
  lambda.0.25se.index <- which(model.cv$lambda == lambda.val.quar.se)
  model.cv$index[4] <- lambda.0.25se.index
  model.cv$index <- c("1se" = lam.1se.index,
                      ".5se" = lambda.0.5se.index,
                      ".25se" = lambda.0.25se.index,
                      "min" = lam.min.index)
  # Which rules have non-zero coef
  nz.coef <- model.cv$nzero[model.cv$index] == 0
  
  if (sum(nz.coef) == 4) {
    stop("All lambda rules without non-zero coefficients. (**)")
  }
  
  
  # while loop to select simplest model if has all zero coefs.
  # simplest model is best because within 1se of min. but make sure not all coefs are zero
  # Preference is lambda.1se. then .5 se. then min
  i <- switch (lambda.sel,
               "1se" = 1,
               "0.5se" = 2,
               "0.25se" = 3,
               "min" = 4
  )
  while (i <= 4) {
    if (nz.coef[i]) {
      i <- i + 1
      next
    }
    model.cv$lambda.value <- model.cv$lambda[model.cv$index[i]]
    model.cv$lambdaRule <- names(model.cv$index)[i]
    break
  }
  
  model.cv$lambdaIndex <- which(model.cv$lambda == model.cv$lambda.value)
  
  bb <- model.cv$glmnet.fit$beta[, model.cv$lambdaIndex]
  model.cv$betas <- bb[bb != 0]
  model.cv$alpha <- 0
  if (predictions) {
    rr <- predict.CVmaxnet(
      object = model.cv,
      newdata = data.frame(data[p == 0, , drop = FALSE]),
      type = "exponent",
      clamp = FALSE,
      offset = offset[p == 0, drop = FALSE])
    raw <- rr/sum(rr)
    model.cv$entropy <- -sum(raw * log(raw))
    model.cv$alpha <- -log(sum(rr))
  } else {
    model.cv$entropy <- NA
    model.cv$alpha <- NA
  }
  model.cv$penalty.factor <- reg
  model.cv$featuremins <- apply(mm, 2, min)
  model.cv$featuremaxs <- apply(mm, 2, max)
  vv <- (sapply(data, class) != "factor")
  model.cv$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model.cv$varmax <- apply(data[, vv, drop = FALSE], 2, max)
  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv], function(n) {
    which.max(table(data[p == 1, n, drop = FALSE]))
  })
  names(majorities) <- names(data)[!vv]
  model.cv$samplemeans <- c(means, majorities)
  model.cv$levels <- lapply(data, levels)
  return(model.cv)
}

CVmaxnet.name <- "CVmaxnet"

# CVmaxnet.fun <- CVmaxnet::CVmaxnet
CVmaxnet.fun <- CVmaxnet

CVmaxnet.errors <- function(occs, envs, bg, tune.args, partitions, algorithm, 
                          partition.settings, other.settings, 
                          categoricals, doClamp, clamp.directions) {
  if(!("rm" %in% names(tune.args)) | !("fc" %in% names(tune.args))) {
    stop("Maxent settings must include 'rm' (regularization multiplier) and 'fc' (feature class) settings. See ?tune.args for details.")
  }else{
    if(!is.numeric(tune.args[["rm"]])) {
      stop("Please input numeric values for 'rm' settings for CVmaxnet.")
    }
    all.fc <- unlist(sapply(1:5, function(x) apply(combn(c("L","Q","H","P","T"), x), 2, function(y) paste(y, collapse = ""))))
    if(any(!tune.args[["fc"]] %in% all.fc)) {
      stop("Please input accepted values for 'fc' settings for CVmaxnet.")
    }
  }
  if(any(tune.args$rm <= 0)) {
    stop("Please input a positive value for 'rm' settings for CVmaxnet.")
  }
}

CVmaxnet.msgs <- function(tune.args, other.settings) {
  # msg <- paste0("CVmaxnet from CVmaxnet package v", packageVersion('CVmaxnet'))
  msg <- "CVmaxnet"
  return(msg)
}

CVmaxnet.args <- function(occs.z, bg.z, tune.tbl.i, other.settings) {
  out <- list()
  out$data <- rbind(occs.z, bg.z)
  out$p <- c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z)))
  out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.tbl.i$fc))
  out$regmult <- tune.tbl.i$rm
  # some models fail to converge if this parameter is not set to TRUE
  # usually the case with sparse datasets
  out$addsamplestobackground <- TRUE
  out <- c(out, other.settings$other.args)
  return(out)
}

CVmaxnet.predict <- function(mod, envs, other.settings) {
  # function to generate a prediction Raster* when raster data is specified as envs,
  # and a prediction data frame when a data frame is specified as envs
  if(inherits(envs, "BasicRaster") == TRUE) {
    envs.n <- raster::nlayers(envs)
    envs.pts <- raster::getValues(envs) %>% as.data.frame()
    mxnet.p <- predict(mod, envs.pts, type = other.settings$pred.type, 
                       clamp = other.settings$doClamp,  other.settings$other.args)
    envs.pts[as.numeric(row.names(mxnet.p)), "pred"] <- mxnet.p
    pred <- raster::rasterFromXYZ(cbind(raster::coordinates(envs), envs.pts$pred), 
                                  res=raster::res(envs), crs = raster::crs(envs)) 
  }else{
    # otherwise, envs is data frame, so return data frame of predicted values
    pred <- predict(mod, envs, type = other.settings$pred.type, na.rm = TRUE, 
                    clamp = other.settings$doClamp, other.settings$other.args) %>% as.numeric()
  }
  return(pred)
}

CVmaxnet.ncoefs <- function(mod) {
  length(mod$betas)
}

# no existing method in model object for variable importance
CVmaxnet.variable.importance <- function(mod) {
  NULL
}

#' @title ENMdetails CVmaxnet
#' @description This is the ENMdetails implementation for CVmaxnet.
#' @export
enm.CVmaxnet <- ENMdetails(name = CVmaxnet.name, fun = CVmaxnet.fun, errors = CVmaxnet.errors,
                         msgs = CVmaxnet.msgs, args = CVmaxnet.args,
                         predict = CVmaxnet.predict, ncoefs = CVmaxnet.ncoefs, variable.importance = CVmaxnet.variable.importance)
