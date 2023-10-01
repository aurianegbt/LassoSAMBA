cpt.crit <- function (cov.model,criterion, method.ll, weight = NULL, pen.coef = NULL) 
{
  ofv <- Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]][["OFV"]]
  ind.model <- Rsmlx:::mlx.getIndividualParameterModel()
  cov.model <- do.call(rbind, cov.model) * 
    1
  if (ncol(cov.model) > 0) {
    cov.info <- Rsmlx:::mlx.getCovariateInformation()
    cat.name <- cov.info$name[grep("categorical", cov.info$type)]
    if (length(cat.name) > 0) {
      cat.cov <- cov.info$covariate[cat.name]
      cat.cov[] <- lapply(cat.cov, factor)
      cat.nb <- unlist(lapply(cat.cov, nlevels)) - 1
      for (j in cat.name) cov.model[, j] <- cov.model[, 
                                                      j] * cat.nb[j]
    }
    i1 <- names(which(ind.model$variability$id))
    i0 <- names(which(!ind.model$variability$id))
    pen.covariate <- sum((cov.model * weight$covariate)[i1, 
    ]) * pen.coef[1] + sum((cov.model * weight$covariate)[i0, 
    ]) * pen.coef[2]
  }
  else {
    pen.covariate <- 0
  }
  cB <- ind.model$correlationBlocks$id
  cor.model <- weight$correlation * 0
  for (k in seq_along(cB)) cor.model[cB[[k]], cB[[k]]] <- 1
  pen.correlation <- sum((lower.tri(cor.model) * cor.model) * 
                           weight$correlation) * pen.coef[1]
  v <- ind.model$variability$id
  pen.variance <- sum(v * weight$variance[names(v)]) * pen.coef[1]
  pen.pop <- length(v) * pen.coef[2]
  error.model <- Rsmlx:::mlx.getContinuousObservationModel()$errorModel
  pen.error <- 0
  for (k in seq_along(error.model)) pen.error <- pen.error + 
    pen.coef[2 + k] * (1 + grepl("combined", error.model[[k]]))
  cr <- ofv + pen.pop + pen.covariate + pen.correlation + 
    pen.variance + pen.error
  return(cr)
}
