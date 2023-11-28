graphsCompTime <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files",project,"/H1.all.R"))

  # Color
  colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR
  # Covariates number arguments
  covariateSizeCar = paste0(covariateSize," covariates")

  # Data to use
  computationStatsCov <- computationStats[computationStats$Method %in% buildMethod
                                          & computationStats$ProjectNumber==covariateSizeCar,]
  ## TIME GRAPHS
  # Value to Display
  valueDisplayTime = data.frame()
  for(t in c("cov","corcov")){
    for(meth in buildMethod){
      aux = computationStatsCov[computationStatsCov$TypeOfSim==t & computationStatsCov$Method==meth,]
      valueDisplayTime = rbind(valueDisplayTime,data.frame(Method=meth,
                                                           TypeOfSim = t,
                                                           Mean = median(aux$time),
                                                           textMean = paste0(round(median(aux$time)), " s"),
                                                           textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                           StandardDeviation = sd(aux$time)))
    }
  }
  valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))


  # Comparative plot
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim %in% c("cov","corcov"),],aes(color=Method,fill=Method,x=factor(Method,levels=buildMethod),y=time))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime, mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates")))+
    #stat_boxplot(geom = "errorbar") +
    xlab("Method used")+
    ylab("Computation time (s)")+
    ggtitle("Computation Time Comparison",
            subtitle=stringr::str_wrap(subtitle,80))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")


  # Save comparative plot
    if(PNG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
             height = 1700, width =  1000+500*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
    }

    if(JPEG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
             height = 1700, width =  1000+500*length(buildMethod), units = "px",device=grDevices::jpeg)
    }

# Cov graphs
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim =="cov",],aes(color=Method,fill=Method,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit")),y=time))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim=="cov",], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    #stat_boxplot(geom = "errorbar") +
    xlab("Method used")+
    ylab("Computation time (s)")+
    ggtitle("Computation Time Comparison",
            subtitle=stringr::str_wrap(paste0(subtitle,"With uncorrelated covariates,"),60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")

  if(PNG){
    ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
            height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
            height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # Corcov graphs
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim =="corcov",],aes(color=Method,fill=Method,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit")),y=time))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim=="corcov",], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    #stat_boxplot(geom = "errorbar") +
    xlab("Method used")+
    ylab("Computation time (s)")+
    ggtitle("Computation Time Comparison",
            subtitle=stringr::str_wrap(paste0(subtitle,"With correlated covariates,"),60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")

  if(PNG){
    ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
            height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
            height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }


  ## ITERATION GRAPHS
  # Value to display
  valueDisplayTime = data.frame()
  for(t in c("cov","corcov")){
    for(meth in buildMethod){
      aux = computationStatsCov[computationStatsCov$TypeOfSim==t & computationStatsCov$Method==meth,]
      valueDisplayTime = rbind(valueDisplayTime,data.frame(Method=meth,
                                                           TypeOfSim = t,
                                                           Mean = median(aux$iteration),
                                                           textMean = paste0(median(aux$iteration)),
                                                           textStandardDeviation = paste0("sd = ", round(sd(aux$iteration),digits=2)),
                                                           StandardDeviation = sd(aux$iteration)))
    }
  }

  valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))

  # Comparative graphs
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim %in% c("cov","corcov"),],aes(color=Method,fill=Method,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit")),y=iteration))+
    geom_violin(lwd=0.5,alpha=0.6)+
    facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates")))+
    geom_text(data=valueDisplayTime, mapping=aes(label=textMean,y=Mean,x=coor,color=Method),vjust=-0.2,fontface="bold",size=6)+
    xlab("Method used")+
    ylab("Iteration Count")+
    ggtitle("Number of Iterations Comparison",
            subtitle=stringr::str_wrap(subtitle,80))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")


  if(PNG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = 1700, width =  1000+500*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }


  if(JPEG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 1700, width =  1000+500*length(buildMethod), units = "px",device=grDevices::jpeg)
  }


  # Cov graphs
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim == "cov",],aes(color=Method,fill=Method,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit")),y=iteration))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim=="cov",], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    xlab("Method used")+
    ylab("Iteration Count")+
    ggtitle("Number of Iterations Comparison",
            subtitle=stringr::str_wrap(paste0(subtitle,"With uncorrelated covariates,"),60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")


  if(PNG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }


  if(JPEG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # Corcov graphs
  ggplot(computationStatsCov[computationStatsCov$TypeOfSim == "corcov",],aes(color=Method,fill=Method,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit")),y=iteration))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim=="corcov",], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    xlab("Method used")+
    ylab("Iteration Count")+
    ggtitle("Number of Iterations Comparison",
            subtitle=stringr::str_wrap(paste0(subtitle,"With correlated covariates,"),60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds"))+
    scale_fill_manual(values=colpas)+
    scale_color_manual(values=colFonce)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")


  if(PNG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
}
