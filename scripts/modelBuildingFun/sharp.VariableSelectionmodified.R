PenalisedRegression <- sharp::PenalisedRegression

sharp.VariableSelection <- function (xdata,
                                     ydata = NULL,
                                     ymeans = NULL,
                                     Lambda = NULL,
                                     pi_list = seq(0.01,0.99, by = 0.01),
                                     K = 100, tau = 0.5, seed = 1,
                                     n_cat = NULL,
                                     family = "gaussian",
                                     implementation = PenalisedRegression, 
                                     resampling = "subsampling",
                                     cpss = FALSE, PFER_method = "MB",
                                     PFER_thr = Inf, FDP_thr = Inf,
                                     Lambda_cardinal = 100, group_x = NULL,
                                     group_penalisation = FALSE,
                                     optimisation = c("grid_search","nloptr"),
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
    if (as.character(substitute(implementation)) %in% c("CART")) {
      Lambda <- seq(1, min(nrow(xdata)/2, 100))
    }
  }
  sharp:::CheckParamRegression(Lambda = Lambda, pi_list = pi_list, 
                       K = K, tau = tau, seed = seed, n_cat = n_cat, family = family, 
                       implementation = implementation, resampling = resampling, 
                       PFER_method = PFER_method, PFER_thr = PFER_thr, FDP_thr = FDP_thr, 
                       Lambda_cardinal = Lambda_cardinal, verbose = verbose)
  
  aux.ydata <- ydata
  
  rownames(xdata) <- paste0("obs",1:nrow(xdata))
  rownames(ymeans) <- paste0("obs",1:nrow(ymeans))
  ydatak = list() 
  for(k in 1:K){
    yaux <- ydata[ydata[,1]==k,2,drop=FALSE]
    rownames(yaux) <- paste0("obs",1:nrow(yaux))
    ydatak <- append(ydatak, list(yaux))
    
  }
  
  xdata.aux <- xdata 
  for(k in 1:K){
    xdata <- xdata.aux 
    ydata <- ydatak[[k]]
    sharp:::CheckDataRegression(xdata = xdata, ydata = ydata, family = family, 
                                verbose = verbose)
    ydatak[[k]] <- ydata
  }
  xdata <- xdata.aux 
  sharp:::CheckDataRegression(xdata = xdata, ydata = ymeans, family = family, 
                              verbose = verbose)
  ymeans <- ydata
  ydata <- aux.ydata
  
  if (group_penalisation) {
    if (is.null(group_x)) {
      stop("Please provide argument 'group_x' for group penalisation. Argument 'group_x' should be a vector with the number of variables in each group.")
    }
  }
  if (is.null(Lambda)) {
    Lambda <- sharp::LambdaGridRegression(xdata = xdata, ydata = ymeans, 
                                   tau = tau, seed = seed, family = family, resampling = resampling, 
                                   Lambda_cardinal = Lambda_cardinal, check_input = FALSE, 
                                   ...)
  }
  optimisation <- match.arg(optimisation)
  extra_args <- list(...)
  
  if (n_cores > 1) {
    if (optimisation != "grid_search") {
      message("Using grid search to allow for parallelisation.")
    }
    future::plan(future::multisession, workers = n_cores)
    
    mypar <- future.apply::future_lapply(X = seq_len(n_cores), 
                                         future.seed = TRUE, FUN = function(k) { 
                                           return(sharp.SerialRegression(xdata = xdata, ydata = ydatak, ymeans=ymeans,
                                                                   Lambda = Lambda,
                                                                   list_rep = split(1:K, cut(1:K, breaks = n_cores, labels = FALSE))[[k]],pi_list = pi_list, K = length(split(1:K, cut(1:K, breaks = n_cores, labels = FALSE))[[k]]), 
                                                                   tau = tau, seed = as.numeric(paste0(seed, 
                                                                                                       k)), n_cat = n_cat, family = family, implementation = implementation, 
                                                                   resampling = resampling, cpss = cpss, PFER_method = PFER_method, 
                                                                   PFER_thr = PFER_thr, FDP_thr = FDP_thr, group_x = group_x, 
                                                                   group_penalisation = group_penalisation, output_data = output_data, 
                                                                   verbose = FALSE, ...))
                                         })
    
    future::plan(future::sequential)
    out <- mypar[[1]]
    if (n_cores > 1) {
      for (i in 2:length(mypar)) {
        out <- do.call(sharp:::Combine, list(stability1 = out, 
                                     stability2 = mypar[[i]]))
      }
    }
  }
  else {
    if (optimisation == "grid_search") {
      out <- sharp.SerialRegression(xdata = xdata, ydata = ydatak,ymeans=ymeans, 
                              Lambda = Lambda, list_rep = 1:K , pi_list = pi_list, K=K , tau = tau, 
                              seed = seed, n_cat = n_cat, family = family, 
                              implementation = implementation, resampling = resampling, 
                              cpss = cpss, PFER_method = PFER_method, PFER_thr = PFER_thr, 
                              FDP_thr = FDP_thr, group_x = group_x, group_penalisation = group_penalisation, 
                              output_data = output_data, verbose = verbose, 
                              ...)
    }
    else {
      eval_f <- function(x, env) {
        out_nloptr <- sharp.SerialRegression(xdata = xdata, 
                                       ydata = ydatak,ymeans=ymeans, Lambda = rbind(x), pi_list = pi_list, list_rep = 1:K, 
                                       K = K, tau = tau, seed = seed, n_cat = n_cat, 
                                       family = family, implementation = implementation, 
                                       resampling = resampling, cpss = cpss, PFER_method = PFER_method, 
                                       PFER_thr = PFER_thr, FDP_thr = FDP_thr, group_x = group_x, 
                                       group_penalisation = group_penalisation, output_data = output_data, 
                                       verbose = FALSE, ...)
        if (any(!is.na(out_nloptr$S))) {
          score <- max(out_nloptr$S, na.rm = TRUE)
        }
        else {
          score <- -Inf
        }
        out <- get("out", envir = env)
        out <- sharp:::Concatenate(out_nloptr, out)
        assign("out", out, envir = env)
        return(-score)
      }
      opts <- list(algorithm = "NLOPT_GN_DIRECT_L", xtol_abs = 0.1, 
                   ftol_abs = 0.1, print_level = 0, maxeval = Lambda_cardinal)
      if ("opts" %in% names(extra_args)) {
        for (opts_id in 1:length(extra_args[["opts"]])) {
          opts[[names(extra_args[["opts"]])[opts_id]]] <- extra_args[["opts"]][[opts_id]]
        }
      }
      nloptr_env <- new.env(parent = emptyenv())
      assign("out", NULL, envir = nloptr_env)
      nloptr_results <- nloptr::nloptr(x0 = apply(Lambda, 
                                                  2, max, na.rm = TRUE), eval_f = eval_f, opts = opts, 
                                       lb = apply(Lambda, 2, min, na.rm = TRUE), ub = apply(Lambda, 
                                                                                            2, max, na.rm = TRUE), env = nloptr_env)
      out <- get("out", envir = nloptr_env)
      out <- sharp:::Concatenate(out, order_output = TRUE)
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
  return(out)
}
