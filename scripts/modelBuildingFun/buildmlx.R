buildmlx <- function(project=NULL,
                     final.project=NULL,
                     model="all",
                     prior=NULL,
                     weight=NULL,
                     coef.w1=0.5,
                     paramToUse="all",
                     covToTest="all",
                     covToTransform="none",
                     center.covariate=FALSE,
                     criterion="BICc",
                     linearization=FALSE,
                     ll=T,
                     test=T,
                     direction=NULL,
                     steps=1000,
                     n.full=10,
                     max.iter=20,
                     explor.iter=2,
                     fError.min=1e-3,
                     seq.cov=FALSE,
                     seq.cov.iter=0,
                     seq.corr=TRUE,
                     p.max=if(buildMethod=="lasso"){1}else{0.1},
                     p.min=c(0.075, 0.05, 0.1),
                     print=TRUE,
                     nb.model=1,
                     nfolds=5,
                     alpha=1,
                     nSS=1000, 
                     buildMethod="lasso",
                     FDR_thr=0.2
                     )
{
  
  doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()))
  ##########################################
  
  ptm <- proc.time()
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"
  op.original <- options()
  op.new <- options()
  op.new$lixoft_notificationOptions$warnings <- 1
  options(op.new)
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data <- resMonolix <- NULL
  pi <- 4 * atan(1)
  if (!is.null(project)) {
    r <- Rsmlx:::prcheck(project, f = "build", paramToUse = paramToUse,
                         model = model)
    if (r$demo)
      return(r$res)
    project <- r$project
  } else {
    project <- Rsmlx:::mlx.getProjectSettings()$project
  }
  method.ll <- iop.ll <- pen.coef <- NULL
  r <- Rsmlx:::buildmlx.check(project, final.project, model, paramToUse,
                              covToTest, covToTransform, center.covariate, criterion,
                              linearization, ll, test, direction, steps, max.iter,
                              explor.iter, seq.cov, seq.cov.iter, seq.corr, p.max,
                              p.min, print, nb.model, prior, weight, n.full)
  if (!is.null(r$change))
    return(list(change = F))
  for (j in 1:length(r)) eval(parse(text = paste0(names(r)[j], "= r[[j]]")))
  r <- Rsmlx:::def.variable(weight = weight, prior = prior, criterion = criterion)
  for (j in 1:length(r)) eval(parse(text = paste0(names(r)[j], "= r[[j]]")))
  is.weight <- weight$is.weight
  is.prior <- NULL
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1",
                   final.project)
  if (dir.exists(final.dir))
    unlink(final.dir, recursive = TRUE)
  project.dir <- Rsmlx:::mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  buildmlx.dir <- file.path(Rsmlx:::mlx.getProjectSettings()$directory,
                            "buildmlx")
  Sys.sleep(0.1)
  if (dir.exists(buildmlx.dir))
    unlink(buildmlx.dir, recursive = TRUE)
  Sys.sleep(0.1)
  dir.create(buildmlx.dir)
  summary.file = file.path(buildmlx.dir, "summary.txt")
  Sys.sleep(0.1)
  if (!dir.exists(final.dir))
    dir.create(final.dir, recursive = T)
  to.cat <- paste0("\n", dashed.line, "\nBuilding:\n")
  if (model$covariate)
    to.cat <- c(to.cat, "  -  The covariate model\n")
  if (model$correlation)
    to.cat <- c(to.cat, "  -  The correlation model\n")
  if (model$residualError)
    to.cat <- c(to.cat, "  -  The residual error model\n")
  to.cat <- c(to.cat, "\n")
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  print.line <- F
  launched.tasks <- Rsmlx:::mlx.getLaunchedTasks()
  p.ini <- Rsmlx:::mlx.getPopulationParameterInformation()
  rownames(p.ini) <- p.ini$name
  ind.omega <- grep("omega_", p.ini[["name"]])
  omega <- p.ini$name[ind.omega]
  omega.ini <- p.ini[ind.omega, ]
  error.model <- Rsmlx:::mlx.getContinuousObservationModel()$errorModel
  obs.dist <- Rsmlx:::mlx.getContinuousObservationModel()$distribution
  covariate.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
  cov.ini <- names(covariate.model[[1]])
  correlation.model <- lapply(Rsmlx:::mlx.getIndividualParameterModel()$correlationBlocks$id,
                              sort)
  if (length(correlation.model) == 0)
    correlation.model <- NULL
  error.model.ini <- error.model
  covariate.model.ini <- covariate.model
  correlation.model.ini <- correlation.model
  to.cat <- ("- - - Initialization - - -\n")
  if (!print.line)
    to.cat <- paste0(plain.line, to.cat)
  Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- Rsmlx:::formatCovariateModel(covariate.model, cov.ini)
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    to.print <- ifelse(!is.null(correlation.model), correlation.model,
                       "NULL")
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- Rsmlx:::formatErrorModel(error.model)
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (!launched.tasks[["populationParameterEstimation"]]) {
    to.cat <- "\nEstimation of the population parameters using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    Rsmlx:::mlx.runPopulationParameterEstimation()
  }
  if (!launched.tasks[["conditionalDistributionSampling"]]) {
    to.cat <- "Sampling of the conditional distribution using the initial model ... \n"
    Rsmlx:::print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    Rsmlx:::mlx.runConditionalDistributionSampling()
  }

  gi <- Rsmlx:::mlx.getSimulatedIndividualParameters()
  gi <- dplyr::filter(gi,rep == gi$rep[nrow(gi)])
  gi = gi[,-which(colnames(gi)=="rep")]
  lin.ll <- method.ll=="linearization"
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]]))  {
      if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
        Rsmlx:::mlx.runConditionalModeEstimation()

      to.cat <- "Estimation of the log-likelihood of the initial model ... \n"
      Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)

      Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
    }

    ll.ini <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
    ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.ini, is.weight, is.prior)
    list.criterion <- ll.ini

    to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
    to.print <- round(ll,2)
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
    #    Rsmlx:::print_result(print, summary.file, to.cat=plain.line, to.print=NULL)

  }

  Rsmlx:::mlx.saveProject(final.project)

  #-----------------------------------------
  if (!is.null(covToTransform)) {
    covariate <- Rsmlx:::mlx.getCovariateInformation()$covariate
    for (cov.name in covToTransform) {
      covkj <- covariate[[cov.name]]
      lckj <- paste0("logt",toupper(substr(cov.name, 1, 1)), substr(cov.name, 2, nchar(cov.name)))
      if (!(lckj %in% Rsmlx:::mlx.getCovariateInformation()$name)) {
        tr.str <- paste0(lckj,' = "log(',cov.name,"/",signif(mean(covkj),digits=2),')"')
        trs <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",tr.str,")")
        eval(parse(text=trs))
        #        colnames(weight$covariate) <- gsub(cov.name, lckj, colnames(weight$covariate))
        foo <- colnames(weight$covariate)
        weight$covariate <- cbind(weight$covariate, weight$covariate[,cov.name])
        colnames(weight$covariate) <- c(foo, lckj)
      }
      covToTest <- c(covToTest, lckj)

    }
    covToTest <- unique(setdiff(covToTest, covToTransform))
    covToTransform <- NULL
    Rsmlx:::mlx.saveProject(final.project)

    if (!Rsmlx:::mlx.getLaunchedTasks()[["populationParameterEstimation"]]) {
      to.cat <- "\nEstimation of the population parameters using the transformed covariates ... \n"
      Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
      Rsmlx:::mlx.runPopulationParameterEstimation()
    }
    if (!Rsmlx:::mlx.getLaunchedTasks()[["conditionalDistributionSampling"]]) {
      to.cat <- "Sampling of the conditional distribution using the the transformed covariates ... \n"
      Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
      Rsmlx:::mlx.runConditionalDistributionSampling()
    }
  }

  #--------------

  stop.test <- FALSE
  corr.test <- FALSE
  iter <- 0

  cov.names0 <- cov.names <- NULL

  if(identical(covToTest,"all")){
    covFix = NULL
  }else{
    covFix <- setdiff(Rsmlx:::mlx.getCovariateInformation()$name, covToTest)
  }

  if (iop.ll) {
    ll.prev <- Inf
    ll.new <- ll.ini
  }
  sp0 <- NULL
  cov.test <- NULL
  e <- NULL
  while (!stop.test) {
    iter <- iter + 1

    if (iter==1){
      weight.covariate <- weight$covariate*coef.w1
      weight.correlation <- weight$correlation*coef.w1
    } else {
      weight.covariate <- weight$covariate
      weight.correlation <- weight$correlation
    }

    to.cat <- paste0(plain.line,"- - - Iteration ",iter," - - -\n")
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
    obs.dist0 <- obs.dist
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll

    if (model$residualError) {
      res.error <- Rsmlx:::errorModelSelection(pen.coef=pen.coef[-c(1, 2)], nb.model=nb.model, f.min=fError.min)
      if (nb.model==1)
        error.model <- res.error
      else {
        error.model <- Rsmlx:::mlx.getContinuousObservationModel()$errorModel
        for (k in (1:length(error.model)))
          error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
      }
    }

    pmax.cov <-  p.max
    # eta = e ; pen.coef = pen.coef[1] ; weight = weight.covariate ; p.max = pmax.cov 
    if (model$covariate){
      if(iter>1){
        lambda.grid = NULL
      }
      ##################################################
      res.covariate <- covariateModelSelection(buildMethod = buildMethod,
                                               nfolds = nfolds,
                                               alpha = alpha,
                                               nSS=nSS,
                                               covFix = covFix,
                                               pen.coef=pen.coef[1],
                                               weight=weight.covariate,
                                               n.full = n.full,
                                               nb.model = nb.model,
                                               direction=direction,
                                               paramToUse=paramToUse,
                                               eta=e,
                                               p.max=pmax.cov,
                                               steps=steps,
                                               sp0=sp0,
                                               iter=iter,
                                               correlation.model=correlation.model,
                                               covariate.model=covariate.model,
                                               criterion = criterion,
                                               FDR_thr=FDR_thr)
      
      ##################################
      res.covariate$res <- Rsmlx:::sortCov(res.covariate$res, cov.ini)
      if (iter>explor.iter) sp0 <- res.covariate$sp
      covToTransform <- setdiff(covToTransform, res.covariate$tr0)
      covariate.model <- res.covariate$model
      e <- res.covariate$residuals
      if (nb.model==1){
        cov.select <- rownames(res.covariate$res)
      }else{
        cov.select <- rownames(res.covariate$res[[1]])}

      cov.names <- lapply(covariate.model[cov.select], function(x) {sort(names(which(x)))})
      cov.names0 <- lapply(covariate.model0[cov.select], function(x) {sort(names(which(x)))})

      cov.test <-Rsmlx:::covariate.test(cov.test, covToTest, covToTransform, paramToUse)

    }else {
      e <- Rsmlx:::mlx.getSimulatedRandomEffects()
    }
    if (model$correlation & !corr.test) {
      if (isTRUE(all.equal(cov.names0,cov.names))) # & isTRUE(all.equal(error.model0,error.model)))
        corr.test <- TRUE
      if (!seq.cov & iter>seq.cov.iter)
        corr.test <- TRUE

      if (corr.test & (seq.cov==T | seq.cov.iter>0)) {
        to.cat <- "Start building correlation model too ... \n"
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
      }
    }
    if (model$correlation & corr.test) {
      pen.corr <- ifelse(iter <= 1, pen.coef[1], pen.coef[1])
      res.correlation <- Rsmlx:::correlationModelSelection(e0=e, pen.coef=pen.corr, nb.model=nb.model,
                                                           corr0=correlation.model0, seqmod=seq.corr, weight=weight.correlation)
      if (nb.model==1)
        correlation.model <- res.correlation
      else
        correlation.model <- res.correlation[[1]]$block
      if (length(correlation.model)==0)
        correlation.model <- NULL
    } else {
      res.correlation <- lapply(Rsmlx:::mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
    }
    if (length(res.correlation)==0)
      res.correlation <- NULL
    #-------------------------------

    if (max.iter>0 | nb.model>1) {
      eq.cov <- isTRUE(all.equal(cov.names0,cov.names)) & (pmax.cov == p.max)
      eq.err <-  isTRUE(all.equal(error.model0,error.model))
      eq.corr <- isTRUE(all.equal(correlation.model0,correlation.model))
      eq.dist <- isTRUE(all.equal(obs.dist0,obs.dist))
      if (!model$correlation | corr.test) {
        if ( eq.cov & eq.err & eq.dist & eq.corr )
          stop.test <- TRUE
        if ( model$covariate & eq.cov & eq.dist & eq.corr )
          stop.test <- TRUE
      }
      if (stop.test) {
        to.cat <- "\nNo difference between two successive iterations\n"
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
      }
      if (!stop.test | nb.model>1) {
        if (model$covariate) {
          to.cat <- "\nCovariate model:\n"
          to.print <- res.covariate$res
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        }
        if (model$correlation) {
          to.cat <- "\nCorrelation model:\n"
          if (!is.null(res.correlation))
            to.print <- res.correlation
          else
            to.print <- "NULL"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        }
        if (model$residualError) {
          to.cat <- "\nResidual error model:\n"
          to.print <- res.error
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        }

      }
    }

    if (!stop.test) {
      p.est <- Rsmlx:::mlx.getEstimatedPopulationParameters()
      Rsmlx:::mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = T)
      p.ini <- Rsmlx:::mlx.getPopulationParameterInformation()
      rownames(p.ini) <- p.ini$name
      i.omega <- which(grepl("omega_", p.ini$name) & (!identical(p.ini$method,"FIXED")))
      p.ini$initialValue[i.omega] <- p.est[p.ini$name[i.omega]]*3
      #          p.ini[omega,] <- omega.ini
      jcor <- grep("corr_",p.ini$name)
      if (length(jcor)>0)  p.ini <- p.ini[-jcor,]
      Rsmlx:::mlx.setPopulationParameterInformation(p.ini)

      if (model$residualError) {
        emodel <- error.model
        odist <- Rsmlx:::mlx.getContinuousObservationModel()$distribution
        for (k in (1:length(emodel))) {
          if (identical(emodel[[k]],"exponential")) {
            emodel[[k]] <- "constant"
            odist[[k]] <- "lognormal"
          } else {
            odist[[k]] <- "normal"
          }
        }
        Rsmlx:::mlx.setErrorModel(emodel)
        Rsmlx:::mlx.setObservationDistribution(odist)
      }

      if (model$covariate) {
        if (length(res.covariate$add.covariate) >0) {
          for (k in 1:length(res.covariate$add.covariate))
            eval(parse(text=res.covariate$add.covariate[[k]]))
        }
        Rsmlx:::mlx.setCovariateModel (covariate.model)
      }

      if (model$correlation & corr.test)
        Rsmlx:::mlx.setCorrelationBlocks(correlation.model)

      #-------------------------------

      if (max.iter>0) {
        if (iop.ll) {
          if (ll.new > ll.prev) {
            g <- Rsmlx:::mlx.getGeneralSettings()
            g$autochains <- FALSE
            g$nbchains <- g$nbchains+1
            Rsmlx:::mlx.setGeneralSettings(g)
          }
          ll.prev <- ll.new
        }
        g= Rsmlx:::mlx.getConditionalDistributionSamplingSettings()
        g$nbminiterations <- max(100, g$nbminiterations)

        #    if (!is.null(g$enablemaxiterations)) {
        # g$enablemaxiterations <- T
        # g$nbmaxiterations <- 500
        #     g$nbminiterations <- 150
        #     g$ratio=0.03
        # }

        Rsmlx:::mlx.setConditionalDistributionSamplingSettings(g)

        # g <- getPopulationParameterEstimationSettings()
        # g$nbexploratoryiterations <-  400
        # g$nbsmoothingiterations <- 200
        # g$smoothingautostop <- g$exploratoryautostop <- F
        # g$simulatedannealing <- T
        # setPopulationParameterEstimationSettings(g)

      }

      Rsmlx:::mlx.saveProject(final.project)
      if (dir.exists(final.dir))
        unlink(final.dir, recursive=TRUE)
      Rsmlx:::mlx.loadProject(final.project)

      if (max.iter>0) {

        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)

        #        Rsmlx:::mlx.runPopulationParameterEstimation(parameter=gi)
        Rsmlx:::mlx.runPopulationParameterEstimation()
        to.cat <- "Sampling from the conditional distribution... \n"
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)

        Rsmlx:::mlx.runConditionalDistributionSampling()
        gi <- Rsmlx:::mlx.getSimulatedIndividualParameters()
        gi <- dplyr::filter(gi,rep==gi$rep[nrow(gi)])
        gi = gi[,-which(colnames(gi)=="rep")]
        if (iop.ll) {
          to.cat <- "Estimation of the log-likelihood... \n"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
            Rsmlx:::mlx.runConditionalModeEstimation()
          Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
          ll.new <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
          list.criterion <- c(list.criterion, ll.new)
        }

        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        Rsmlx:::mlx.saveProject(buildmlx.project.iter)

        if (iop.ll) {
          if (stop.test)
            ll <- ll0
          else {
            ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight)
          }

          to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
          to.print <- round(ll,2)
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        }
        if (iter >= max.iter) {
          stop.test <- TRUE
          to.cat <- "Maximum number of iterations reached\n"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
        }
      }
    }
    if (max.iter==0)
      stop.test <- TRUE
    Rsmlx:::mlx.saveProject(final.project)
  }
  change.error.model <- NULL
  if (iop.ll) {
    ll.min <- min(list.criterion)
    iter.opt <- which.min(list.criterion)
    if (iter.opt==1)
      buildmlx.project.iter <- project
    else {
      buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter.opt-1,".mlxtran"))
    }
    Rsmlx:::mlx.loadProject(buildmlx.project.iter)
    Rsmlx:::mlx.saveProject(final.project)
  } else {
    iter.opt <- iter
  }

  if (test) {

    if (!iop.ll) {
      g <- as.list(Rsmlx:::mlx.getLaunchedTasks())$logLikelihoodEstimation
      if (!linearization & !("importanceSampling" %in% g))
        Rsmlx:::mlx.runLogLikelihoodEstimation(linearization=F)
      if (linearization & !("linearization" %in% g))
        Rsmlx:::mlx.runLogLikelihoodEstimation(linearization=T)
      ll.min <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
    }

    to.cat <- paste0(plain.line,"- - - Further tests - - -\n")
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
    if (model$covariate & sum(Rsmlx:::mlx.getIndividualParameterModel()$variability$id)>1) {
      g0 <- Rsmlx:::mlx.getIndividualParameterModel()
      covariate <- random.effect <- p.ttest <- p.lrt <- in.model <- p.value <- NULL
      r.test <-Rsmlx:::covariate.test(cov.test, covToTest, covToTransform, paramToUse)
      r.test <- dplyr::select(dplyr::filter(r.test,!in.model),-in.model)
      if (is.weight) {
        w.cov <- weight$covariate[cbind(r.test[['parameter']], r.test[['covariate']])]
        r.test <- dplyr::mutate(r.test,p.value = Rsmlx:::p.weight(p.value, w.cov, pen.coef[1]))
      }
      r.cov0 <- res.covariate$r.cov0
      for (j in 1:nrow(r.test)) {
        pj <- r.test$parameter[j]
        if (!is.null(r.cov0[[pj]])) {
          if (r.test$covariate[j] %in% r.cov0[[pj]])
            r.test$p.value[j] <- 1
        }
      }
      i.min <- which(as.numeric(r.test$p.value) < p.min[1])
      if (length(i.min)>0) {

        to.cat <- paste0(plain.short,"Add parameters/covariates relationships:\n")
        to.print <- r.test[i.min,]
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        g <- Rsmlx:::mlx.getIndividualParameterModel()
        stop.test <- F
        for (i in i.min) {
          param.i <- r.test$parameter[i]
          cov.i <- r.test$covariate[i]
          g$covariateModel[[param.i]][cov.i] <- T
        }
        Rsmlx:::mlx.setIndividualParameterModel(g)
        iter <- iter+1
        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        Rsmlx:::mlx.saveProject(buildmlx.project.iter)
        Rsmlx:::mlx.runPopulationParameterEstimation()
      } else {
        stop.test <- F
      }
      if (any(Rsmlx:::mlx.getObservationInformation()$type != "continuous"))
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=F)
      else
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=T)
      #Rsmlx:::mlx.runStandardErrorEstimation(linearization = lin.ll)

      # if (is.weight | is.prior) {
      #   w.cov <- weight$covariate[cbind(r.test[['parameter']], r.test[['covariate']])]
      #   r.test <- r.test %>%  mutate(p.value =Rsmlx:::p.weight(p.value, w.cov, pen.coef[1]))
      # }
      g <- Rsmlx:::mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])

      r.test <- Rsmlx:::mlx.getTests()$wald
      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0.99999
      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            g$covariateModel[[np]][nc] <- F
            ipc <- grep(paste0("beta_",np,"_",nc), r.test$parameter)
            pv[ipc] <- Rsmlx:::p.weight(pv[ipc], weight$covariate[np, nc], pen.coef[1])
            if (min(pv[ipc]) < p.min[2])
              g$covariateModel[[np]][nc] <- T
            else
              list.ipc <- c(list.ipc, ipc)
          }
        }
      }

      if (identical(g$covariateModel, g0$covariateModel))
        stop.test <- T

      if (length(list.ipc) >0) {
        to.cat <- paste0(plain.short,"Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (dplyr::rename(dplyr::select(r.test,-c(method, statistics)),coefficient=parameter))[list.ipc,]
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
      }

      g1 <- g
      if (!stop.test) {
        if (length(list.ipc) >0) {
          Rsmlx:::mlx.setIndividualParameterModel(g)
          iter <- iter+1
          to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
          Rsmlx:::mlx.saveProject(buildmlx.project.iter)
          Rsmlx:::mlx.runPopulationParameterEstimation()
        }
        if (lin.ll) {
          if(!launched.tasks[["conditionalModeEstimation"]])
            Rsmlx:::mlx.runConditionalModeEstimation()
        } else {
          to.cat <- "Sampling from the conditional distribution... \n"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          Rsmlx:::mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
        Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
        ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
        to.print <- round(ll,2)
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        if (ll.new < ll.min) {
          ll.min <- ll.new
          Rsmlx:::mlx.saveProject(final.project)
        } else {
          Rsmlx:::mlx.loadProject(final.project)
        }
      }

      # ---   Wald tests
      g <- as.list(Rsmlx:::mlx.getLaunchedTasks())$standardErrorEstimation
      if (!linearization & !("stochasticApproximation" %in% g))
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=F)
      if (!("linearization" %in% g))
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=T)

      r.test <- Rsmlx:::mlx.getTests()$wald
      g <- Rsmlx:::mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])

      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0

      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            ipc <- grep(paste0("beta_",np,"_",nc), r.test$parameter)
            pv[ipc] <- Rsmlx:::p.weight(pv[ipc], weight$covariate[np, nc], pen.coef[1])
            if (max(pv[ipc]) > p.min[2]) {
              g$covariateModel[[np]][nc] <- F
              list.ipc <- c(list.ipc, ipc)
            }
          }
        }
      }

      if (length(list.ipc) > 0 & !identical(g$covariateModel, g0$covariateModel) & !identical(g$covariateModel, g1$covariateModel))  {
        to.cat <- paste0(plain.short,"Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (dplyr::rename(dplyr::select(r.test,-c(method, statistics)),coefficient=parameter))[list.ipc,]
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)

        Rsmlx:::mlx.setIndividualParameterModel(g)
        iter <- iter+1
        to.cat <- paste0("\nRun scenario for model ",iter," ... \nEstimation of the population parameters... \n")
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
        buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
        Rsmlx:::mlx.saveProject(buildmlx.project.iter)
        Rsmlx:::mlx.runPopulationParameterEstimation()
        if (lin.ll) {
          if(!launched.tasks[["conditionalModeEstimation"]])
            Rsmlx:::mlx.runConditionalModeEstimation()
        } else {
          to.cat <- "Sampling from the conditional distribution... \n"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          Rsmlx:::mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
        Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
        ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
        to.print <- round(ll,2)
        Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
        if (ll.new < ll.min) {
          ll.min <- ll.new
          Rsmlx:::mlx.saveProject(final.project)
        }
      }

      # ---  end Wald tests

    }
    if (model$correlation) {
      test.cor <- T
      cor.block0 <- Rsmlx:::mlx.getIndividualParameterModel()$correlationBlocks$id
      while (test.cor) {
        Rsmlx:::mlx.loadProject(final.project)
        p.cortest <- NULL
        if (!Rsmlx:::mlx.getLaunchedTasks()$conditionalDistributionSampling)
          Rsmlx:::mlx.runConditionalDistributionSampling()
        r.test <- dplyr::rename( dplyr::filter(Rsmlx:::correlationTest()$p.value,!in.model),p.value=p.cortest)
        param1 <- gsub("eta_","",r.test$randomEffect.1)
        param2 <- gsub("eta_","",r.test$randomEffect.2)
        w.cor <- weight$correlation[cbind(param1, param2)]+weight$correlation[cbind(param2, param1)]
        r.test <- dplyr::mutate(r.test ,p.value =Rsmlx:::p.weight(p.value, w.cor, pen.coef[1]))

        i.min <- which(as.numeric(r.test$p.value) < p.min[3])
        g <- Rsmlx:::mlx.getIndividualParameterModel()
        param.list <- unlist(g$correlationBlocks$id)
        if (length(i.min)>0) {
          i.min <- i.min[which.min(r.test$p.value[i.min])]
          param1 <- gsub("eta_","",r.test$randomEffect.1[i.min])
          param2 <- gsub("eta_","",r.test$randomEffect.2[i.min])
          test.cor <- F
          if ( !(param1 %in% param.list) & !(param2 %in% param.list) ) {
            l.block <- length(g$correlationBlocks$id)+1
            g$correlationBlocks$id[[l.block]] <- c(param1, param2)
            test.cor <- T
          }
          if ( !(param1 %in% param.list) & (param2 %in% param.list) ) {
            l.block <- grep(param2, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]], param1)
            test.cor <- T
          }
          if ( (param1 %in% param.list) & !(param2 %in% param.list) ) {
            l.block <- grep(param1, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]], param2)
            test.cor <- T
          }
          if (test.cor) {

            to.cat <- paste0(plain.short, "Add correlation:\n")
            to.print <- r.test[i.min,]
            Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
            to.print <- g$correlationBlocks$id
            Rsmlx:::print_result(print, summary.file, to.cat=NULL, to.print=to.print)
            gi <- Rsmlx:::mlx.getPopulationParameterInformation()
            gi$initialValue[which(gi$name==paste0("omega_",param1))] <- 3*gi$initialValue[which(gi$name==paste0("omega_",param1))]
            gi$initialValue[which(gi$name==paste0("omega_",param2))] <- 3*gi$initialValue[which(gi$name==paste0("omega_",param2))]
            Rsmlx:::mlx.setPopulationParameterInformation(gi)
            Rsmlx:::mlx.setIndividualParameterModel(g)

            ###
            iter <- iter+1
            buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
            Rsmlx:::mlx.saveProject(buildmlx.project.iter)
            to.cat <- paste0("\nRun scenario for model ",iter,"  ... \nEstimation of the population parameters... \n")
            Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
            Rsmlx:::mlx.runPopulationParameterEstimation()
            if (lin.ll) {
              if(!launched.tasks[["conditionalModeEstimation"]])
                Rsmlx:::mlx.runConditionalModeEstimation()
            } else {
              to.cat <- "Sampling from the conditional distribution... \n"
              Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
              Rsmlx:::mlx.runConditionalDistributionSampling()
            }
            to.cat <- "Estimation of the log-likelihood... \n"
            Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
            Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
            ll.new <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
            ll.disp <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)

            to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
            to.print <- round(ll.disp,2)
            Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)

            if (ll.new < ll.min) {
              ll.min <- ll.new
              Rsmlx:::mlx.saveProject(final.project)
            } else {
              test.cor <- F
            }

          } else {
            test.cor <- F
          }
        } else {
          test.cor <- F
        }
      }

      Rsmlx:::mlx.loadProject(final.project)
      p <- Rsmlx:::mlx.getEstimatedPopulationParameters()
      if (any(Rsmlx:::mlx.getObservationInformation()$type != "continuous")) {
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=F)
        se <- Rsmlx:::mlx.getEstimatedStandardErrors()$stochasticApproximation$se
        names(se) <- Rsmlx:::mlx.getEstimatedStandardErrors()$stochasticApproximation$parameter
      } else {
        Rsmlx:::mlx.runStandardErrorEstimation(linearization=T)
        se <- Rsmlx:::mlx.getEstimatedStandardErrors()$linearization$se
        names(se) <- Rsmlx:::mlx.getEstimatedStandardErrors()$linearization$parameter
      }
      z <- as.numeric(p)/as.numeric(se[names(p)])
      names(z) <- names(p)
      pv <- pnorm(-abs(z))*2
      pv.corr <- pv[grep("corr_", names(pv))]
      if (length(which(pv.corr>p.min[2])) > 0) {
        Rsmlx:::mlx.saveProject(buildmlx.project.iter)
        ind <- Rsmlx:::mlx.getEstimatedIndividualParameters()$saem
        pv.block <- strsplit(gsub("corr_", "", names(pv.corr)),"_")
        ind.mod <- Rsmlx:::mlx.getIndividualParameterModel()
        cb <- ind.mod$correlationBlocks$id
        cl <- unlist(cb)
        pv.tot <- rep(0,length(cl))
        names(pv.tot) <- cl
        for (k in 1: length(pv.corr))
          pv.tot[pv.block[[k]]] <- pv.tot[pv.block[[k]]] + log(pv.corr[k])
        param0 <- names(which.max(pv.tot))
        cb <- lapply(cb, function(x) setdiff(x, param0))
        i.cb <- which(lapply(cb, length)==1)
        if (length(i.cb)>0)
          cb[[i.cb]] <- NULL
        if (!identical(cor.block0, cb)) {
          ind.mod$correlationBlocks$id <- cb
          to.cat <- paste0(plain.short,"Test correlation model:\n")
          to.print <- cb
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)

          Rsmlx:::mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = TRUE)
          Rsmlx:::mlx.setIndividualParameterModel(ind.mod)
          iter <- iter+1
          buildmlx.project.iter <- file.path(buildmlx.dir,paste0("iteration",iter,".mlxtran"))
          Rsmlx:::mlx.saveProject(buildmlx.project.iter)
          to.cat <- paste0("Run scenario for model ",iter,"  ... \nEstimation of the population parameters... \n")
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          Rsmlx:::mlx.runPopulationParameterEstimation(parameters=ind)
          if (lin.ll) {
            if(!launched.tasks[["conditionalModeEstimation"]])
              Rsmlx:::mlx.runConditionalModeEstimation()
          } else {
            to.cat <- "Sampling from the conditional distribution... \n"
            Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
            Rsmlx:::mlx.runConditionalDistributionSampling()
          }
          to.cat <- "Estimation of the log-likelihood... \n"
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)
          Rsmlx:::mlx.runLogLikelihoodEstimation(linearization = lin.ll)
          ll.new <- Rsmlx:::compute.criterion(criterion, method.ll, weight, pen.coef)
          ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], criterion, ll.new, is.weight, is.prior)
          to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
          to.print <- round(ll,2)
          Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
          if (ll.new < ll.min) {
            ll.min <- ll.new
            Rsmlx:::mlx.saveProject(final.project)
          }
        }
      }
    }

  }

  if (model$covariate & nb.model>1)
    res.covariate$res <- sortCov(res.covariate$res[[1]], cov.ini)

  Rsmlx:::mlx.loadProject(final.project)

  if (model$covariate)
    covariate.model.print <- Rsmlx:::formatCovariateModel(Rsmlx:::mlx.getIndividualParameterModel()$covariateModel)
  if (model$correlation) {
    correlation.model.print <- lapply(Rsmlx:::mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
    if (length(correlation.model.print)==0)
      correlation.model.print <- NULL
  }
  if (model$residualError)
    error.model.print <- Rsmlx:::formatErrorModel(Rsmlx:::mlx.getContinuousObservationModel()$errorModel)
  if (iop.ll) {
    ind.model <- Rsmlx:::mlx.getIndividualParameterModel()
    cov.model <- do.call(rbind, ind.model$covariateModel) *  1
    print(cov.model)
    ll.final <- Rsmlx:::compute.criterion(criterion, method.ll, 
                                  weight, pen.coef)
    ll <- Rsmlx:::formatLL(Rsmlx:::mlx.getEstimatedLogLikelihood()[[method.ll]], 
                   criterion, ll.final, is.weight, is.prior)
  }

  to.cat <- paste0(plain.line,"\nFinal statistical model:\n")
  Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)

  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- covariate.model.print
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    if (!is.null(correlation.model.print))
      to.print <- correlation.model.print
    else
      to.print <- "NULL"
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- error.model.print
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
  }
  if (iop.ll & max.iter>0) {
    to.cat <- paste0("\nEstimated criteria (",method.ll,"):\n")
    to.print <- round(ll,2)
    Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=to.print)
  }

  test.del <- FALSE
  if (model$covariate & center.covariate) {
    foo <- lapply(res.covariate$model,function(x) {which(x)})
    cov.model <- unique(unlist(lapply(foo,function(x) {names(x)})))
    cov.type <- Rsmlx:::mlx.getCovariateInformation()$type[cov.model]
    cov.cont <- names(cov.type[cov.type=="continuous"])
    for (ck in cov.cont) {
      cck <- paste0("c",ck)
      covk <- Rsmlx:::mlx.getCovariateInformation()$covariate[[ck]]
      tr.str <- paste0(cck,' = "',ck,"-",signif(mean(covk),digits=2),'"')
      tr.str <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",tr.str,")")
      eval(parse(text=tr.str))
      g= Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
      cg <- lapply(g, function(x) {foo <- x[ck]; x[ck]<-x[cck]; x[cck]=foo; return(x)})
      Rsmlx:::mlx.setCovariateModel (cg)
      test.del <- TRUE
      covariate.model <- cg
    }
  }
  if (test.del & dir.exists(final.dir)) {
    unlink(final.dir, recursive=TRUE)
  }

  g= Rsmlx:::mlx.getScenario()
  g$tasks[[2]]=TRUE
  Rsmlx:::mlx.setScenario(g)
  Rsmlx:::mlx.saveProject(final.project)
  Rsmlx:::mlx.loadProject(final.project)

  error.model <- Rsmlx:::mlx.getContinuousObservationModel()$errorModel
  covariate.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
  correlation.model <- lapply(Rsmlx:::mlx.getIndividualParameterModel()$correlationBlocks$id, sort)
  if (length(correlation.model)==0)
    correlation.model <- NULL

  dt <- proc.time() - ptm
  res <- list(project=final.project, niter=iter.opt, time=dt["elapsed"])
  if (model$covariate)
    res <- c(res, list(covariate.model=covariate.model))
  if (model$correlation)
    res <- c(res, list(correlation.model=correlation.model))
  if (model$residualError)
    res <- c(res, list(error.model=error.model))

  ttime = round(dt["elapsed"], digits=1)
  to.cat <- paste0("\ntotal time: ", round(dt["elapsed"], digits=1),"s\n", plain.line)
  Rsmlx:::print_result(print, summary.file, to.cat=to.cat, to.print=NULL)


  res$change <- !(identical(error.model,error.model.ini) &
                    identical(covariate.model,covariate.model.ini) &
                    identical(correlation.model, correlation.model.ini))
  res$change.error.model <- change.error.model
  res$weight <- weight
  res$covToTest <- covToTest
  options(op.original)
  
  #################################
  return(list(res=res,time=ttime,iter=iter))
  #################################
}
