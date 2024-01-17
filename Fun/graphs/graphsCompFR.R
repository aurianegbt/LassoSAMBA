graphsCompFR <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Color
  colFonce = colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR

  # Data to use
  errorStatsCov <- errorStats[errorStats$Method %in% buildMethod,]
  errorStatsParCov <-  suppressMessages(errorStatsPar %>%
    group_by(Model,ProjectNumber,TypeOfSim,Method) %>%
    summarise(
      FP = sum(FP),
      TP = sum(TP), 
      FN = sum(FN),
      TN = sum(TN)
      ) %>% 
      mutate(FDR = (FP/max(FP+TP,1)),.after = "TP") %>%
      mutate(FNR = (FN/max(FN+TN,1)),.after="TN") %>%
      as.data.frame())


  # FDR
  # Value to Display
  valueDisplayFDR = data.frame()
  for(t in c("cov","corcov")){
    for(meth in buildMethod){
      aux = errorStatsParCov[errorStatsParCov$TypeOfSim==t & errorStatsParCov$Method==meth,]
      valueDisplayFDR = rbind(valueDisplayFDR,data.frame(Method=meth,
                                                         TypeOfSim = t,
                                                         Mean = mean(aux$FDR),
                                                         Median = median(aux$FDR),
                                                         text = paste0("Mean : ",round(mean(aux$FDR),digits=2),"\nMedian : ",round(median(aux$FDR),digits=2)),
                                                         textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                         StandardDeviation = sd(aux$time)))
    }
  }

  valueDisplayFDR <- cbind(valueDisplayFDR,coor=1:length(buildMethod))

  # Comparative plot
  ggplot(errorStatsCov[errorStatsCov$TypeOfSim %in% c("cov","corcov"),],aes(color=Method,fill=Method,y=FDR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassoSSnoCov0","elasticnetSSnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayFDR, mapping=aes(label=text,x=coor,color=Method),y=0.10,vjust=-0.2,fontface="bold",size=3)+
    facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates"))) +
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=stringr::str_wrap(subtitle,90))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
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
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = 1700, width = 1000+500*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 1700, width = 1000+500*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # Cov plot
  ggplot(errorStatsCov[errorStatsCov$TypeOfSim =="cov",],aes(color=Method,fill=Method,y=FDR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassoSSnoCov0","elasticnetSSnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayFDR[valueDisplayFDR$TypeOfSim=="cov",], mapping=aes(label=text,x=coor,color=Method),y=0.10,vjust=-0.2,fontface="bold",size=3)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
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
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = 1700, width =  800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # Corcov plot
  ggplot(errorStatsCov[errorStatsCov$TypeOfSim =="corcov",],aes(color=Method,fill=Method,y=FDR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassoSSnoCov0","elasticnetSSnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayFDR[valueDisplayFDR$TypeOfSim=="corcov",], mapping=aes(label=text,x=coor,color=Method),y=0.10,vjust=-0.2,fontface="bold",size=3)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
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
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorFDR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # # FNR
  # # Value to Display
  # valueDisplayFNR = data.frame()
  # for(t in c("cov","corcov")){
  #   for(meth in buildMethod){
  #     aux = errorStatsParCov[errorStatsParCov$TypeOfSim==t & errorStatsParCov$Method==meth,]
  #     valueDisplayFNR = rbind(valueDisplayFNR,data.frame(Method=meth,
  #                                                        TypeOfSim = t,
  #                                                        Mean = mean(aux$FNR),
  #                                                        Median = median(aux$FNR),
  #                                                        text = paste0("Mean : ",round(mean(aux$FNR),digits=2),"\nMedian : ",round(median(aux$FNR),digits=5)),
  #                                                        textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
  #                                                        StandardDeviation = sd(aux$time)))
  #   }
  # }
  # 
  # valueDisplayFNR <- cbind(valueDisplayFNR,coor=1:length(buildMethod))
  # 
  # # Comparative plot
  # ggplot(errorStatsCov[errorStatsCov$TypeOfSim %in% c("cov","corcov"),],aes(color=Method,fill=Method,y=FNR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"))))+
  #   geom_violin(lwd=0.5,alpha=0.6)+
  #   geom_text(data=valueDisplayFNR, mapping=aes(label=text,x=coor,color=Method),y=0.0075,vjust=-0.2,fontface="bold",size=3)+
  #   facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates",cov2="With uncorrelated covariates",corcov2="With correlated covariates"))) +
  #   guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
  #   xlab("Method used")+
  #   ggtitle("False Negative Rate Distribution",
  #           subtitle=stringr::str_wrap(subtitle,90))+
  #   scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
  #   scale_fill_manual(values=colpas)+
  #   scale_color_manual(values=colFonce)+
  #   theme(axis.text.x = element_text(size = 10))+
  #   theme(axis.text.y = element_text(size = 8))+
  #   theme(axis.title = element_text(size=12))+
  #   theme(strip.text = element_text(size = 12))+
  #   theme(plot.subtitle = element_text(size=12))+
  #   theme(plot.title = element_text(size=16,color="#ee6c4d"))+
  #   theme(legend.position = "none")
  # 
  # if(PNG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
  #          height = 1700, width =  1000+500*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  # }
  # 
  # if(JPEG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
  #          height = 1700, width =  1000+500*length(buildMethod), units = "px",device=grDevices::jpeg)
  # }
  # 
  # # Cov plot
  # ggplot(errorStatsCov[errorStatsCov$TypeOfSim =="cov",],aes(color=Method,fill=Method,y=FNR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"))))+
  #   geom_violin(lwd=0.5,alpha=0.6)+
  #   geom_text(data=valueDisplayFNR[valueDisplayFNR$TypeOfSim=="cov",], mapping=aes(label=text,x=coor,color=Method),y=0.0075,vjust=-0.2,fontface="bold",size=3)+
  #   guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
  #   xlab("Method used")+
  #   ggtitle("False Negative Rate Distribution",
  #           subtitle=stringr::str_wrap(subtitle,60))+
  #   scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
  #   scale_fill_manual(values=colpas)+
  #   scale_color_manual(values=colFonce)+
  #   theme(axis.text.x = element_text(size = 10))+
  #   theme(axis.text.y = element_text(size = 8))+
  #   theme(axis.title = element_text(size=12))+
  #   theme(strip.text = element_text(size = 12))+
  #   theme(plot.subtitle = element_text(size=12))+
  #   theme(plot.title = element_text(size=16,color="#ee6c4d"))+
  #   theme(legend.position = "none")
  # 
  # if(PNG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
  #          height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  # }
  # 
  # if(JPEG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
  #          height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  # }
  # 
  # # Corcov plot
  # ggplot(errorStatsCov[errorStatsCov$TypeOfSim =="corcov",],aes(color=Method,fill=Method,y=FNR,x=factor(Method,levels=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"))))+
  #   geom_violin(lwd=0.5,alpha=0.6)+
  #   geom_text(data=valueDisplayFNR[valueDisplayFNR$TypeOfSim=="corcov",], mapping=aes(label=text,x=coor,color=Method),y=0.0075,vjust=-0.2,fontface="bold",size=3)+
  #   guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
  #   xlab("Method used")+
  #   ggtitle("False Negative Rate Distribution",
  #           subtitle=stringr::str_wrap(subtitle,60))+
  #   scale_x_discrete(labels=c(reg="stepAIC",lassoSS="Lasso",elasticnetSS="Elastic Net",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test",lassoSSREP="Lasso with\ns.s. on replicates",elasticnetSSREP="Elastic Net with\ns.s. on replicates",lassoSSCritREP="Lasso with mult.\nthresholds and s.s. on rep.",elasticnetSSCritREP="Elastic Net with mult.\nthresholds and s.s. on rep."))+
  #   scale_fill_manual(values=colpas)+
  #   scale_color_manual(values=colFonce)+
  #   theme(axis.text.x = element_text(size = 10))+
  #   theme(axis.text.y = element_text(size = 8))+
  #   theme(axis.title = element_text(size=12))+
  #   theme(strip.text = element_text(size = 12))+
  #   theme(plot.subtitle = element_text(size=12))+
  #   theme(plot.title = element_text(size=16,color="#ee6c4d"))+
  #   theme(legend.position = "none")
  # 
  # if(PNG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
  #          height = 1700, width =   800+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  # }
  # 
  # if(JPEG){
  #   ggsave(paste0(Folder,"/ErrorFNR",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
  #          height = 1700, width =   800+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  # }
}
