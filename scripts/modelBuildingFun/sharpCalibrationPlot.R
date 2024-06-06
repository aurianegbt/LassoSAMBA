sharpCalibrationPlot.init <- function(project=NULL,pathToSave=NULL,title.param=NULL,nSS=1000,alpha=1,nfolds=5,exclude=NULL,FDP_thr=0.05,JPEG=F,PNG=T){
  
  library(ggplot2)

  doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()))
  
  suppressMessages({
    if (!is.null(project)){
      project <- Rsmlx:::prcheck(project)$project
    }else{
      project <- Rsmlx:::mlx.getProjectSettings()$project
    }
    
    lixoftConnectors::loadProject(project)
  })
  
  if(is.null(pathToSave)){
    pathToSave = paste0(stringr::str_sub(lixoftConnectors::getProjectSettings()$directory,end = -7),"/sharpCalibrationPlot_")
  }
  
  if(!lixoftConnectors::getLaunchedTasks()$populationParameterEstimation){
    lixoftConnectors::runPopulationParameterEstimation()
  }
  if(!lixoftConnectors::getLaunchedTasks()$conditionalDistributionSampling){
    lixoftConnectors::runConditionalDistributionSampling()
  }
  
  lixoftConnectors::saveProject()
  
  sp.df <- Rsmlx:::mlx.getSimulatedIndividualParameters()
  if (is.null(sp.df$rep))
    sp.df$rep <- 1
  nrep <- max(sp.df$rep)
  ind.dist <- Rsmlx:::mlx.getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  n.param <- length(param.names)
  cov.info <- Rsmlx:::mlx.getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  tcov.names <- NULL
  covariates <- cov.info$covariate
  cov.cat <- cov.names[cov.types == "categorical"]
  covariates[cov.cat] <- lapply(covariates[cov.cat], as.factor)
  indvar <- Rsmlx:::mlx.getIndividualParameterModel()$variability$id
  
  cov.model <- Rsmlx:::mlx.getIndividualParameterModel()$covariateModel
  eps <- 1e-15
  
  # Y transformation
  Y=sp.df[,c("rep","id",names(indvar)[which(indvar==TRUE)])]

  for(j in 1:n.param){
    dj <- ind.dist[j]
    nj <- names(dj)
    if (indvar[j]) {
      if (tolower(dj) == "lognormal") {
        Y[,nj] <- log(Y[,nj] + eps)
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("log.", nj)
      }else if (tolower(dj) == "logitnormal") {
        Y[,nj] <- log((Y[,nj] + eps)/(1 - Y[,nj] + eps))
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("logit.", nj)
      }else if (tolower(dj) == "probitnormal") {
        Y[,nj] <- qnorm(Y[,nj])
        colnames(Y)[which(colnames(Y)==nj)]  <- paste0("probit.", nj)
      }
    }
  }
  
  Sigma=diag(Rsmlx:::mlx.getEstimatedPopulationParameters()[paste0("omega_",param.names[which(indvar)])]**2)
  colnames(Sigma) <- param.names[which(indvar)]
  rownames(Sigma) <- param.names[which(indvar)]
  #######" FAIRE LA SELECTION ICI 
  N = length(unique(Y$id))
  
  Y.mat = sapply(Y[,-c(1,2)],function(x){rowMeans(matrix(x,nrow=N))}) #1 : rep 2 : id
  X.mat  = covariates[,setdiff(colnames(covariates),"id")]
  
  library(foreach)
  
  
  foreach(p = names(indvar)[which(indvar)],.packages=c("sharp","ggplot2","grDevices","gghighlight"),.export = "CalibrationPlot") %dopar% { 
    source("~/Travail/00_Theme.R")
    plot = CalibrationPlot(Y.mat[,stringr::str_detect(colnames(Y.mat),p),drop=F],
                    X.mat,
                    Sigma[p,p],
                    nfolds,alpha,nSS,p.name=p,title.param=title.param[[p]],
                    n_cores = floor(parallel::detectCores()/sum(which(indvar))),
                    FDP_thr = FDP_thr)
    if(PNG){
      ggsave(plot=plot$plot1, filename =paste0(pathToSave,p,".png"),
             height=1500,width=3000,units = "px",bg='transparent',device=grDevices::png,create.dir = TRUE)
    }
    if(JPEG){
      ggsave(plot=plot$plot1, filename =paste0(pathToSave,p,".jpeg"),
             height=1500,width=3000,units = "px",device=grDevices::jpeg,create.dir = TRUE)
    }
    
    # if(PNG){
    #   ggsave(plot=plot$plot2, filename =paste0(pathToSave,p,"15.png"),
    #          height=1500,width=3000,units = "px",bg='transparent',device=grDevices::png)
    # }
    # if(JPEG){
    #   ggsave(plot=plot$plot2, filename =paste0(pathToSave,p,"15.jpeg"),
    #          height=1500,width=3000,units = "px",device=grDevices::jpeg)
    # }
  }
}

CalibrationPlot <- function(Y,X,omega,
                            nfolds=5,
                            alpha=1,
                            nSS=1000,
                            p.name=NULL,
                            title.param=NULL,
                            n_cores=NULL,
                            FDP_thr=0.05){
  
  cov.names = colnames(X)
  tparam.names = colnames(Y)
  
  if(!is.matrix(Y)){
    Yaux <- as.matrix(Y)
  }else{Yaux=Y}
  if(!is.matrix(X)){
    Xaux <- as.matrix(X)
  }else{Xaux=X}
  
  Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
  
  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  Ywh <- Yaux %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)
  colnames(Xwh) <- cov.names
  
  VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,nfolds=nfolds,alpha=alpha,K=nSS,n_cores=n_cores,FDP_thr = FDP_thr)
  
  pi_list = VariableSelection.outputs$params$pi_list
  lambda_list = VariableSelection.outputs$Lambda
  
  
  Score = VariableSelection.outputs$S_2d
  df = data.frame()
  for(i in 1:ncol(Score)){
    df <- rbind(df,data.frame(lambda = signif(lambda_list,digits=2),pi = pi_list[i],Score=Score[,i]))
  }
  
  max = setNames(as.list(sharp::Argmax(VariableSelection.outputs)),c("lambda","pi"))
  
  browser()
  
  plot = ggplot(df,aes(x=as.factor(lambda),y=pi,fill=Score)) +
    geom_tile() +
    # geom_segment(x=as.factor(min(df$lambda)),xend=as.factor(max(df$lambda)),y=max$pi,yend=max$pi,linewidth=1,color="black") +
    # geom_segment(y=min(df$pi),yend=max(df$pi),x=as.factor(max$lambda),xend=max$lambda,linewidth=1,color="black") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_gradientn(colours = c("ivory", "navajowhite", "tomato","darkred"),na.value="white") + 
    xlab(latex2exp::TeX("$\\lambda$"))+ylab(latex2exp::TeX("$t_{SS}$")) +
    theme(axis.title=element_text(size=14,face="bold")) +
    theme(axis.text.x=element_text(size=8)) +
    scale_y_continuous(breaks=c(seq(0,0.5,0.05),seq(0.6,1,0.1)))  + 
    ggtitle(label=title.param,
            subtitle=latex2exp::TeX(r"(Heatmap of the stability score S$(t_{SS},\lambda)$)")) + 
    theme(plot.subtitle = element_text(size=10)) + 
    theme(legend.position = "bottom") + 
    theme(legend.text = element_text(face="plain"))  +
    theme(legend.title = element_text(size=10))
  
  plot2 = ggplot(df,aes(x=as.factor(lambda),y=pi,fill=Score)) +
    geom_tile() +
    gghighlight::gghighlight(Score >= quantile(df$Score,0.85,na.rm=T),unhighlighted_params = list(fill="grey")) +
    # geom_segment(x=as.factor(min(df$lambda)),xend=as.factor(max(df$lambda)),y=max$pi,yend=max$pi,linewidth=1,color="black") +
    # geom_segment(y=min(df$pi),yend=max(df$pi),x=as.factor(max$lambda),xend=max$lambda,linewidth=1,color="black") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_gradientn(colours = c("ivory", "navajowhite", "tomato","darkred"),na.value="white") + 
    xlab(latex2exp::TeX("$\\lambda$"))+ylab(latex2exp::TeX("$t_{SS}$")) +
    theme(axis.title=element_text(size=14,face="bold")) +
    theme(axis.text.x=element_text(size=8)) +
    scale_y_continuous(breaks=c(seq(0,0.5,0.05),seq(0.6,1,0.1)))  + 
    ggtitle(label=title.param,
            subtitle=latex2exp::TeX(r"(Heatmap of the stability score S$(t_{SS},\lambda)$ of 15% higher scores.)")) + 
    theme(plot.subtitle = element_text(size=10)) + 
    theme(legend.position = "bottom") + 
    theme(legend.text = element_text(face="plain"))  +
    theme(legend.title = element_text(size=10))
  
  
  return(list(plot1=plot,plot2=plot2))
}
