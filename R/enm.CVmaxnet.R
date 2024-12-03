################################# #
# maxnet ENMdetails object ####
################################# #

CVmaxnet.name <- "CVmaxnet"

CVmaxnet.fun <- CVmaxnet::CVmaxnet

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
  msg <- paste0("CVmaxnet from CVmaxnet package v", packageVersion('CVmaxnet'))
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
  out$verbose <- FALSE
  out$predictions <- FALSE
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
