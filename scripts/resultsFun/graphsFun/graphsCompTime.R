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
                                                         Median = median(aux$time), 
                                                         Q76 = quantile(aux$time,0.76),
                                                         textMean = round(median(aux$time)),
                                                         textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                         StandardDeviation = sd(aux$time)))
  }
  valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))
  
  # graphs
  plot <- ggplot(computationStatsCov,aes(color=Method,fill=Method,x=factor(Method,levels=buildMethod),y=time))+
    geom_boxplot(lwd=0.5,alpha=0.6,width=0.5)+
    geom_text(data=valueDisplayTime, mapping=aes(label=textMean,y=Median,x=coor-0.25),vjust=-0.2,hjust=1,fontface="bold",size=4)+
    xlab("Method used")+
    ylab("Computation time (s)")+
    scale_x_discrete(labels=c(setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))],"lassoFDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))]),
                              stepAIC="stepAIC\nwith stat. test",
                              SAEMVS="SAEMVS"))+
    # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_fill_manual(values=setNames(colpas,buildMethod))+
    scale_color_manual(values=setNames(colFonce,buildMethod))+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(legend.position = "none")
  
  for(k in 1:length(buildMethod)){
    plot <- plot + geom_segment(lwd=1,x=valueDisplayTime[valueDisplayTime$Method==buildMethod[k],"coor"]-0.45,xend=valueDisplayTime[valueDisplayTime$Method==buildMethod[k],"coor"],y=valueDisplayTime[valueDisplayTime$Method==buildMethod[k],"Median"],yend=valueDisplayTime[valueDisplayTime$Method==buildMethod[k],"Median"],color=colFonce[k])
  }
  
  
  if(PNG){
    ggsave(paste0(Folder,"/ComputationTime.png"),
           height = 800, width =   500+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ComputationTime.jpeg"),
           height = 800, width =   500+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  
  ggplot(computationStatsCov,aes(color=Method,fill=Method,x=factor(Method,levels=buildMethod),y=iteration))+
    geom_boxplot(lwd=0.5,alpha=0.6)+
    xlab("Method used")+
    ylab("Iteration")+
    # ggtitle("Computation Time Comparison",
            # subtitle=stringr::str_wrap(subtitle,60))+
    scale_x_discrete(labels=c(stepAIC="stepAIC\nwith stat. test",
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))],"lassoFDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))]),
                              SAEMVS="SAEMVS"))+
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
    ggsave(paste0(Folder,"/Iter.png"),
           height = 800, width =   500+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/Iter.jpeg"),
           height = 800, width =   500+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
}
