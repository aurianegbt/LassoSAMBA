graphsLL <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){
  
  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Color
  colFonce = colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR

  # Data to use
  likelihoodStatsCov <- likelihoodStats[likelihoodStats$Method %in% buildMethod,]

  plot <- 
    ggplot(likelihoodStatsCov,aes(x=Criterion,y=Value,fill=Method))+geom_boxplot()+
    scale_fill_manual(values=setNames(colpas,buildMethod),
                      labels=c(stepAIC="stepAIC\nwith stat. test",
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))],"lassoFDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))]),
                              SAEMVS="SAEMVS"))+
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size=12),
          strip.text = element_text(size = 12),
          plot.subtitle = element_text(size=12),
          legend.position="bottom",
          legend.text = element_text(size=10),
          legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))
    
  
  if(PNG){
    ggsave(paste0(Folder,"/ICcomparison.png"),
           height = 800, width =  500*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ICcomparison.jpeg"),
           height = 800, width =   500*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  return(plot)
}
