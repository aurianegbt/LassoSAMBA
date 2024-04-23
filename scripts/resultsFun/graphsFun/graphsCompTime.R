graphsCompTime <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Color
  colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR
  
  # Data to use
  computationStatsCov <- computationStats[computationStats$Method %in% buildMethod,]
  
  ## TIME GRAPHS
  # Value to Display
  valueDisplayTime = data.frame()
  for(meth in buildMethod){
    aux = computationStatsCov[ computationStatsCov$Method==meth,]
    valueDisplayTime = rbind(valueDisplayTime,data.frame(Method=meth,
                                                         Mean = median(aux$time),
                                                         textMean = paste0(round(median(aux$time)), " s"),
                                                         textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                         StandardDeviation = sd(aux$time)))
  }
  valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))

# graphs
  ggplot(computationStatsCov,aes(color=Method,fill=Method,x=factor(Method,levels=buildMethod),y=time))+
    geom_violin(lwd=0.5,alpha=0.6)+
    geom_text(data=valueDisplayTime, mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6)+
    xlab("Method used")+
    ylab("Computation time (s)")+
    ggtitle("Computation Time Comparison",
            subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(reg="stepAIC\nwhith stat. test",lasso="lasso\nwhithout s.s.",elastinet="elastic net\nwhithout s.s.",lassoSS="Lasso\nwhith stat. test",elasticnetSS="Elastic Net\nwhith stat. test",rlasso="Lasso with\ns.s. on replicates",relasticnet="Elastic Net with\ns.s. on replicates",lassoCrit = "Lasso with\nmultiple thresholds and no s.s.",elasticnetCrit="Elastic Net with\nmultiple thresholds and no s.s.",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",rlassoCrit="Lasso with mult.\nthresholds and s.s. on rep.",relasticnetCrit="Elastic Net with mult.\nthresholds and s.s. on rep.",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC",lassoSSCov0="Lasso",elasticnetSSnoCov0="Elastic Net"))+
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
    ggsave(paste0(Folder,"/ComputationTime.png"),
            height = 1700, width =   1200+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ComputationTime.jpeg"),
            height = 1700, width =   1200+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
}
