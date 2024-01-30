VariableSelection <- function (xdata, ydata = NULL,
                               Lambda = NULL,
                               pi_list = seq(0.01,0.99, by = 0.01),
                               K = 100,
                               tau = 0.5,
                               seed = 1,
                               n_cat = NULL, 
                               family = "gaussian",
                               implementation = PenalisedRegression, 
                               resampling = "subsampling", 
                               cpss = FALSE, 
                               PFER_method = "MB", 
                               PFER_thr = Inf,
                               FDP_thr = Inf, 
                               Lambda_cardinal = 100,
                               group_x = NULL, 
                               group_penalisation = FALSE,
                               n_cores = 1,
                               output_data = FALSE, 
                               verbose = TRUE,
                               beep = NULL, ...) 
{
  if (is.null(Lambda)) {
    if (as.character(substitute(implementation)) %in% c("SparseGroupPLS", 
                                                        "GroupPLS")) {
      Lambda <- seq(1, length(group_x) - 1)
    }
    if (as.character(substitute(implementation)) %in% c("SparsePLS", 
                                                        "SparsePCA")) {
      Lambda <- seq(1, ncol(xdata) - 1)
    }
  }
  sharp:::CheckParamRegression(Lambda = Lambda, pi_list = pi_list, 
                       K = K, tau = tau, seed = seed, n_cat = n_cat, family = family, 
                       implementation = implementation, resampling = resampling, 
                       PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr, 
                       Lambda_cardinal = Lambda_cardinal, verbose = verbose)
  # for(k in 1:length(ydata)){
  #   sharp:::CheckDataRegression(xdata = xdata, ydata = ydata[[k]], family = family, 
  #                               verbose = verbose)
  # }
  if (group_penalisation) {
    if (is.null(group_x)) {
      stop("Please provide argument 'group_x' for group penalisation. Argument 'group_x' should be a vector with the number of variables in each group.")
    }
  }
  if (is.null(Lambda)) {
    Lambda <- sharp:::LambdaGridRegression(xdata = xdata, ydata = ydata[[1]], 
                                   tau = tau, seed = seed, family = family, resampling = resampling, 
                                   Lambda_cardinal = Lambda_cardinal, check_input = FALSE,...)
  }
 mypar <- SerialRegression(xdata = xdata, ydata = ydata, 
                            Lambda = Lambda, pi_list = pi_list, K = ceiling(K/n_cores), 
                            tau = tau, seed = as.numeric(paste0(seed,rpois(1,2))), n_cat = n_cat, 
                            family = family, implementation = implementation, 
                            resampling = resampling, cpss = cpss, PFER_method = PFER_method, 
                            PFER_thr = PFER_thr, FDP_thr = FDP_thr, group_x = group_x, 
                            group_penalisation = group_penalisation, output_data = output_data, 
                            verbose = verbose,...)
  out <- mypar
  if (n_cores > 1) {
    for (i in 2:length(mypar)) {
      out <- do.call(Combine, list(stability1 = out, stability2 = mypar[[i]]))
    }
  }
  if ("methods" %in% names(out)) {
    myimplementation <- as.character(substitute(implementation))
    if (is.function(resampling)) {
      myresampling <- as.character(substitute(resampling))
    }
    else {
      myresampling <- resampling
    }
    out$methods$implementation <- myimplementation
    out$methods$resampling <- myresampling
  }
  class(out) <- "variable_selection"
  if (!is.null(beep)) {
    beepr::beep(sound = beep)
  }
  return(sharp::SelectedVariables(out))
}
