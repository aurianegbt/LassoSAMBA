graphsCompFR <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Color
  colFonce = colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR

  # Data to use
  errorStatsCov <- errorStats[errorStats$Method %in% buildMethod,]
  errorStatsParCov <-  suppressMessages(errorStatsPar %>%
    group_by(Model,Method) %>%
    summarise(
      FP = sum(FP),
      TP = sum(TP), 
      FN = sum(FN),
      TN = sum(TN)
      ) %>% 
      mutate(FDR = (FP/max(FP+TP,1)),.after = "TP") %>%
      mutate(FNR = (FN/max(FN+TN,1)),.after="TN")%>%
      mutate(F1_score = TP/(TP+1/2*(FN+FP)),.after="FNR") %>%
      as.data.frame())
  
  errorStatsParCov <- errorStatsParCov[errorStatsParCov$Method %in% buildMethod,]


  # FDR
  # Value to Display
  valueDisplayFDR = data.frame()
  for(meth in buildMethod){
    aux = errorStatsParCov[errorStatsParCov$Method==meth,]
    valueDisplayFDR = rbind(valueDisplayFDR,data.frame(Method=meth,
                                                       Mean = mean(aux$FDR),
                                                       Median = median(aux$FDR),
                                                       text = paste0("Mean : ",round(mean(aux$FDR),digits=2),"\nMedian : ",round(median(aux$FDR),digits=2)),
                                                       textStandardDeviation = paste0("sd = ", round(sd(aux$FDR),digits=2)),
                                                       StandardDeviation = sd(aux$FDR)))
  }
  
  valueDisplayFDR <- cbind(valueDisplayFDR,coor=1:length(buildMethod))
  
  # Cov plot
  plotFDR = ggplot(errorStatsParCov,aes(color=Method,fill=Method,y=FDR,x=factor(Method,levels=c("reg","lasso","elasticnet","lassoSS","elasticnetSS","rlasso","relasticnet","lassoCrit","elasticnetCrit","lassoSSCrit","elasticnetSSCrit","rlassoCrit","reslaticnetCrit","regnoCov0","lassoSSnoCov0","lassoSSCritnoCov0","elasticnetSSnoCov0",buildMethod[stringr::str_detect(buildMethod,"regPEN")],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp","sharpnoCov0"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayFDR, mapping=aes(label=text,x=coor),color="black",y=0.10,vjust=-0.2,fontface="bold",size=3)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC\nwith stat. test",
                              lasso="lasso\nwithout s.s.",
                              elastinet="elastic net\nwithout s.s.",
                              lassoSS="Lasso\nwith stat. test",
                              elasticnetSS="Elastic Net\nwith stat. test",
                              rlasso="Lasso with\ns.s. on replicates", 
                              relasticnet="Elastic Net with\ns.s. on replicates",   
                              lassoCrit = "Lasso with\nmultiple thresholds and no s.s.",  
                              elasticnetCrit="Elastic Net with\nmultiple thresholds and no s.s.",  
                              lassoSSCritnoCov0="Lasso with\nadapted initialisation",   
                              elasticnetSSCrit="Elastic Net with\nmultiple thresholds", 
                              rlassoCrit="Lasso with mult.\nthresholds and s.s. on rep.",
                              relasticnetCrit="Elastic Net with mult.\nthresholds and s.s. on rep.",   
                              regnoCov0="stepAIC",     
                              lassoSSnoCov0="Lasso", 
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],"sharpnoCov0"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))]),
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
                              
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))]),
                              setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),
                              elasticnetSSnoCov0="Elastic Net",
                              sharp = "Lasso calibrated using sharp\nwith stat. test",
                              sharpnoCov0 = "Lasso calibrated using sharp"))+
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
    ggsave(paste0(Folder,"/ErrorFDR.png"),
           height = 1700, width =  1000+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ErrorFDR.jpeg"),
           height = 1700, width =   1000+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }

  # FNR plot
  # Value to Display
  valueDisplayFNR = data.frame()
  for(meth in buildMethod){
    aux = errorStatsParCov[errorStatsParCov$Method==meth,]
    valueDisplayFNR = rbind(valueDisplayFNR,data.frame(Method=meth,
                                                       Mean = mean(aux$FNR),
                                                       Median = median(aux$FNR),
                                                       text = paste0("Mean : ",round(mean(aux$FNR),digits=2),"\nMedian : ",round(median(aux$FNR),digits=2)),
                                                       textStandardDeviation = paste0("sd = ", round(sd(aux$FNR),digits=2)),
                                                       StandardDeviation = sd(aux$FNR)))
  }
  
  valueDisplayFNR <- cbind(valueDisplayFNR,coor=1:length(buildMethod))
  
  # Cov plot
  plotFNR = ggplot(errorStatsParCov,aes(color=Method,fill=Method,y=FNR,x=factor(Method,levels=c("reg","lasso","elasticnet","lassoSS","elasticnetSS","rlasso","relasticnet","lassoCrit","elasticnetCrit","lassoSSCrit","elasticnetSSCrit","rlassoCrit","reslaticnetCrit","regnoCov0","lassoSSnoCov0","lassoSSCritnoCov0","elasticnetSSnoCov0",buildMethod[stringr::str_detect(buildMethod,"regPEN")],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp","sharpnoCov0"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayFNR, mapping=aes(label=text,x=coor),color="black",y=0.0010,vjust=-0.2,fontface="bold",size=3)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Negative Rate Distribution",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC\nwith stat. test",
                              lasso="lasso\nwithout s.s.",
                              elastinet="elastic net\nwithout s.s.",
                              lassoSS="Lasso\nwith stat. test",
                              elasticnetSS="Elastic Net\nwith stat. test",
                              rlasso="Lasso with\ns.s. on replicates", 
                              relasticnet="Elastic Net with\ns.s. on replicates",   
                              lassoCrit = "Lasso with\nmultiple thresholds and no s.s.",  
                              elasticnetCrit="Elastic Net with\nmultiple thresholds and no s.s.",  
                              lassoSSCritnoCov0="Lasso with\nadapted initialisation",   
                              elasticnetSSCrit="Elastic Net with\nmultiple thresholds", 
                              rlassoCrit="Lasso with mult.\nthresholds and s.s. on rep.",
                              relasticnetCrit="Elastic Net with mult.\nthresholds and s.s. on rep.",   
                              regnoCov0="stepAIC",     
                              lassoSSnoCov0="Lasso", 
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],"sharpnoCov0"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))]),
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
                              
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))]),
                              setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),
                              elasticnetSSnoCov0="Elastic Net",
                              sharp = "Lasso calibrated using sharp\nwith stat. test",
                              sharpnoCov0 = "Lasso calibrated using sharp"))+
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
    ggsave(paste0(Folder,"/ErrorFNR.png"),
           height = 1700, width =  1100+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ErrorFNR.jpeg"),
           height = 1700, width =   1100+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  
  # F1-score plot 
  # Value to Display
  valueDisplayF1_score = data.frame()
  for(meth in buildMethod){
    aux = errorStatsParCov[errorStatsParCov$Method==meth,]
    valueDisplayF1_score = rbind(valueDisplayF1_score,data.frame(Method=meth,
                                                       Mean = mean(aux$F1_score),
                                                       Median = median(aux$F1_score),
                                                       text = paste0("Mean : ",round(mean(aux$F1_score),digits=2),"\nMedian : ",round(median(aux$F1_score),digits=2)),
                                                       textStandardDeviation = paste0("sd = ", round(sd(aux$F1_score),digits=2)),
                                                       StandardDeviation = sd(aux$F1_score)))
  }
  
  valueDisplayF1_score <- cbind(valueDisplayF1_score,coor=1:length(buildMethod))
  
  # Cov plot
  plotF1_score = ggplot(errorStatsParCov,aes(color=Method,fill=Method,y=F1_score,x=factor(Method,levels=c("reg","lasso","elasticnet","lassoSS","elasticnetSS","rlasso","relasticnet","lassoCrit","elasticnetCrit","lassoSSCrit","elasticnetSSCrit","rlassoCrit","reslaticnetCrit","regnoCov0","lassoSSnoCov0","lassoSSCritnoCov0","elasticnetSSnoCov0",buildMethod[stringr::str_detect(buildMethod,"regPEN")],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp","sharpnoCov0"))))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayF1_score, mapping=aes(label=text,x=coor),color="black",y=0.8,vjust=-0.2,fontface="bold",size=3)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("F1-score Distribution",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC\nwith stat. test",
                              lasso="lasso\nwithout s.s.",
                              elastinet="elastic net\nwithout s.s.",
                              lassoSS="Lasso\nwith stat. test",
                              elasticnetSS="Elastic Net\nwith stat. test",
                              rlasso="Lasso with\ns.s. on replicates", 
                              relasticnet="Elastic Net with\ns.s. on replicates",   
                              lassoCrit = "Lasso with\nmultiple thresholds and no s.s.",  
                              elasticnetCrit="Elastic Net with\nmultiple thresholds and no s.s.",  
                              lassoSSCritnoCov0="Lasso with\nadapted initialisation",   
                              elasticnetSSCrit="Elastic Net with\nmultiple thresholds", 
                              rlassoCrit="Lasso with mult.\nthresholds and s.s. on rep.",
                              relasticnetCrit="Elastic Net with mult.\nthresholds and s.s. on rep.",   
                              regnoCov0="stepAIC",     
                              lassoSSnoCov0="Lasso", 
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],"sharpnoCov0"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))]),
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))]),
                              setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),
                              elasticnetSSnoCov0="Elastic Net",
                              sharp = "Lasso calibrated using sharp\nwith stat. test",
                              sharpnoCov0 = "Lasso calibrated using sharp"))+
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
    ggsave(paste0(Folder,"/ErrorF1.png"),
           height = 1700, width =  1100+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ErrorF1.jpeg"),
           height = 1700, width =   1100+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  
  plotFDR <- plotFDR  + ggtitle("",subtitle="")
  plotFNR <- plotFNR  + ggtitle("",subtitle="")
  plotF1_score <- plotF1_score  + ggtitle("",subtitle="")
  
  gp <- 
    annotate_figure(ggarrange(plotFDR,plotFNR,plotF1_score, ncol =3),
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ) %>%
    annotate_figure(
      top=text_grob("Error Statistics Distribution",
                    face="bold",size=18,color="#862B0D")
    )
  
  if(PNG){
    ggsave(paste0(Folder,"/Error.png"),
           height = 1500, width =  3000+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/Error.jpeg"),
           height = 1500, width =   3000+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  
}
