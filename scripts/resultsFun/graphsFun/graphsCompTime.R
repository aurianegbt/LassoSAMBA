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
    scale_x_discrete(labels=c(reg="stepAIC with stat. test",
                              lasso="Lasso without s.s.",
                              elastinet="Elastic net without s.s.",
                              lassoSS="Lasso with stat. test",
                              elasticnetSS="Elastic Net with stat. test",
                              rlasso="Lasso with s.s. on replicates",
                              relasticnet="Elastic Net with s.s. on replicates",
                              lassoCrit = "Lasso with multiple thresholds and no s.s.",
                              elasticnetCrit="Elastic Net with multiple thresholds and no s.s.",
                              lassoSSCrit="Lasso with multiple thresholds",
                              elasticnetSSCrit="Elastic Net with multiple thresholds",
                              rlassoCrit="Lasso with mult. thresholds and s.s. on rep.",
                              relasticnetCrit="Elastic Net with mult. thresholds and s.s. on rep.",
                              setNames(paste0("penalized stepAIC pen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),
                              regnoCov0="stepAIC",
                              lassoSSnoCov0="Lasso",
                              lassoSSCritnoCov0="Lasso",
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],"sharpnoCov0"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))]),
                              setNames(paste0("Lasso\nFDP<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))]),
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
                              elasticnetSSnoCov0="Elastic Net",
                              sharpnoCov0="Lasso calibrated using sharp",
                              sharp="Lasso calibrated using sharp\nwith stat. test"))+
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
