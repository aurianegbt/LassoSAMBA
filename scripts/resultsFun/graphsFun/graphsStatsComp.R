graphsStatsComp <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){
  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))
  
  # Color
  colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR
  
  # Data.frame to use
  resultModelParCov <- resultModelPar[resultModelPar$Method %in% buildMethod, ]
  errorStatsParCov <- errorStatsPar[errorStatsPar$Method %in% buildMethod,]
  errorStatsParCov <-  suppressMessages(errorStatsParCov %>%
    group_by(Model,Method) %>%
    summarise(
      FP = sum(FP),
      TP = sum(TP), 
      FN = sum(FN),
      TN = sum(TN)
    )) %>%
    as.data.frame()

  errorStatsParCov <- errorStatsParCov  %>% 
    mutate(FDR = (FP/sapply(FP+TP,FUN=function(x){max(x,1)})),.after = "TP") %>%
    mutate(FNR = (FN/sapply(FN+TN,FUN=function(x){max(x,1)})),.after="TN") %>%
    mutate(F1_score = TP/(TP+1/2*(FN+FP)),.after="FNR")
  
  value = data.frame()
  for(meth in buildMethod){
    value <- rbind(value,data.frame(Method=meth,
                                    Type=c("Exact","Overselection"),
                                    Value=c(nrow(errorStatsParCov[errorStatsParCov$Method==meth & errorStatsParCov$FP==0 & errorStatsParCov$FN==0,]),nrow(errorStatsParCov[errorStatsParCov$Method==meth  & errorStatsParCov$FN==0,]))))
  }
  
  plot <- ggplot(mapping=aes(x=factor(Method,levels=buildMethod), y=Value,pattern=Type,fill = Method)) +
    geom_col(data=value[value$Type=="Overselection",], aes(group = Method), width = 0.7, position = position_dodge(width = 0.7),color="black")  +
    geom_bar_pattern(data=value[value$Type=="Exact",], position="dodge",width=0.7,color="black",pattern_fill="black",pattern_density=0.05,pattern_spacing=0.025,stat='identity') +
    labs(x="Method", y = "Proportion (in %)") +
    scale_fill_brewer(palette="Set2")+ 
    scale_y_continuous(breaks=seq(0,100,10),limits = c(0,100))+ 
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
                              setNames(paste0("Lasso\nE[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))])[order(as.numeric(stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP")))],
                              setNames(paste0("Lasso\n",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
                              elasticnetSSnoCov0="Elastic Net",
                              sharpnoCov0="Lasso calibrated using sharp",
                              sharp="Lasso calibrated using sharp\nwith stat. test")) +
    scale_pattern_manual(values = c(Exact="stripe","Overselection"="none"))+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    guides(fill = "none")+
    theme(legend.position = "bottom")
  
  if(PNG){
    ggsave(paste0(Folder,"/ComparisonStats.png"),
           height = 1200, width =   500+300*length(buildMethod), units = "px", bg='transparent',device=grDevices::png)
  }
  
  if(JPEG){
    ggsave(paste0(Folder,"/ComparisonStats.jpeg"),
           height = 1200, width =   500+300*length(buildMethod), units = "px",device=grDevices::jpeg)
  }
  
}
