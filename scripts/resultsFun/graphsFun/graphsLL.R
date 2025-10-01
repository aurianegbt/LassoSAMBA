graphsLL <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){
  
  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Color
  colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43","#9C4400","#005E57","#7A003C")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FF7F11","#00B8A9","#D81159")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#FFD2A6","#BFF0E6","#F7B8D2")[1:length(buildMethod)]
  

  # Data to use
  likelihoodStatsCov <- likelihoodStats[likelihoodStats$Method %in% buildMethod,]

  plot <- 
    ggplot(likelihoodStatsCov,aes(x=Criterion,y=Value,fill=Method))+geom_boxplot()+
    scale_fill_manual(values=setNames(colpas,buildMethod),
                      labels=c(stepAIC="step-SAMBA",
                               setNames(paste0("lasso-SAMBA\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))],"lassoFDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"lassoFDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"lassoFDP"))]),
                               setNames(paste0("Elastic Net\nalpha=",ifelse(stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"elasticnet") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"elasticnet"))],"elasticnet")==10,"1",ifelse(stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"elasticnet") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"elasticnet"))],"elasticnet")==0,"0",paste0("0.",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"elasticnet") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"elasticnet"))],"elasticnet"))))),buildMethod[stringr::str_detect(buildMethod,"elasticnet") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"elasticnet"))]),
                               SAEMVS="SAEMVS",lassoBIC="lassoBIC-SAMBA",lassoSSrepFDP10="lassorep-SAMBA\nE[FDR]<10%"))+
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
