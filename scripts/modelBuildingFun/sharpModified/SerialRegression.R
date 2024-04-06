SerialRegression <- function (xdata, ydata = NULL, Lambda, pi_list = seq(0.6, 0.9, 
                                                     by = 0.01), K = 100, tau = 0.5, seed = 1, n_cat = 3, family = "gaussian", 
          implementation = PenalisedRegression, resampling = "subsampling", 
          cpss = FALSE, PFER_method = "MB", PFER_thr = Inf, FDP_thr = Inf, 
          group_x = NULL, group_penalisation = FALSE, output_data = FALSE, 
          verbose = TRUE, ...) 
{
  N <- N_block <- ncol(xdata)
  Xsub <- xdata
  Ysub <- ydata[[1]]
  
  
  mybeta <- sharp::SelectionAlgo(xdata = Xsub, ydata = Ysub, Lambda = Lambda[, 
                                                                      1], group_x = group_x, family = family, implementation = implementation,...)
  Beta <- array(0, dim = c(nrow(mybeta$selected), ncol(mybeta$selected), 
                           length(ydata)))
  rownames(Beta) <- rownames(mybeta$selected)
  colnames(Beta) <- colnames(mybeta$selected)
  if (length(dim(mybeta$beta_full)) == 2) {
    Beta_full <- array(0, dim = c(dim(mybeta$beta_full)[1], 
                                  dim(mybeta$beta_full)[2], length(ydata)), dimnames = list(rownames(mybeta$beta_full), 
                                                                                dimnames(mybeta$beta_full)[[2]], NULL))
  }else {
    if (length(dim(mybeta$beta_full)) == 3) {
      Beta_full <- array(0, dim = c(dim(mybeta$beta_full)[1], 
                                    dim(mybeta$beta_full)[2], length(ydata), dim(mybeta$beta_full)[3]), 
                         dimnames = list(rownames(mybeta$beta_full), 
                                         dimnames(mybeta$beta_full)[[2]], NULL, dimnames(mybeta$beta_full)[[3]]))
    }else {
      stop(paste0("Invalid output from the variable selection function provided in 'implementation'. The output 'beta_full' must be an array with 2 or 3 dimensions."))
    }
  }
  withr::local_seed(seed)
  
  allbeta = foreach(k = 1:length(ydata)) %dopar%{
    Xsub <- xdata
    Ysub <- ydata[[k]]
    mybeta <- sharp::SelectionAlgo(xdata = Xsub, ydata = Ysub, 
                                   Lambda = Lambda[, 1], group_x = group_x, family = family, 
                                   implementation = implementation,...)
    return(mybeta)
  }
  for(k in 1:length(ydata)){
    Beta[rownames(allbeta[[k]]$selected), colnames(allbeta[[k]]$selected), 
         k] <- allbeta[[k]]$selected
    if (length(dim(Beta_full)) == 3) {
      Beta_full[rownames(allbeta[[k]]$beta_full), colnames(allbeta[[k]]$beta_full), 
                k] <- allbeta[[k]]$beta_full
    }
    else {
      Beta_full[rownames(allbeta[[k]]$beta_full), colnames(allbeta[[k]]$beta_full), 
                k, ] <- allbeta[[k]]$beta_full
    }
  }
  
  
  
  bigstab <- matrix(NA, nrow = nrow(Beta), ncol = ncol(Beta))
  colnames(bigstab) <- colnames(Beta)
  rownames(bigstab) <- rownames(Beta)
  for (i in 1:nrow(Beta)) {
    for (j in 1:ncol(Beta)) {
      bigstab[i, j] <- sum(Beta[i, j, ] != 0, na.rm = TRUE)/sum(!is.na(Beta[i, 
                                                                            j, ]))
    }
  }
  if (group_penalisation) {
    metrics <- sharp::StabilityMetrics(selprop = bigstab, pk = NULL, 
                                pi_list = pi_list, K = length(ydata), n_cat = n_cat, Sequential_template = NULL, 
                                graph = FALSE, group = group_x, PFER_method = PFER_method, 
                                PFER_thr_blocks = PFER_thr, FDP_thr_blocks = FDP_thr)
  } else {
    metrics <- sharp::StabilityMetrics(selprop = bigstab, pk = NULL, 
                                pi_list = pi_list, K = length(ydata), n_cat = n_cat, Sequential_template = NULL, 
                                graph = FALSE, PFER_method = PFER_method, PFER_thr_blocks = PFER_thr, 
                                FDP_thr_blocks = FDP_thr)
  }
  Beta <- Beta_full
  myimplementation <- as.character(substitute(implementation, 
                                              env = parent.frame(n = 2)))
  if (is.function(resampling)) {
    myresampling <- as.character(substitute(resampling))
  } else {
    myresampling <- resampling
  }
  out <- list(S = metrics$S, Lambda = Lambda, Q = metrics$Q, 
              Q_s = metrics$Q_s, P = metrics$P, PFER = metrics$PFER, 
              FDP = metrics$FDP, S_2d = metrics$S_2d, PFER_2d = metrics$PFER_2d, 
              FDP_2d = metrics$FDP_2d, selprop = bigstab, Beta = Beta, 
              methods = list(type = "variable_selection", implementation = myimplementation, 
                             family = family, resampling = myresampling, cpss = cpss, 
                             PFER_method = PFER_method), params = list(K = length(ydata), 
                                                                       pi_list = pi_list, tau = tau, n_cat = n_cat, pk = ncol(xdata), 
                                                                       n = nrow(xdata), PFER_thr = PFER_thr, FDP_thr = FDP_thr, 
                                                                       seed = seed))
  if (output_data) {
    out$params <- c(out$params, list(xdata = xdata, ydata = ydata))
  }
  class(out) <- "variable_selection"
  return(out)
}
