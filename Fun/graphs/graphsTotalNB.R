graphsTotalNB <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files",project,"/H1.all.R"))

  # Color & covariates
  gr = "#888888"
  fill.vec = c(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(Reduce(union,H1.all))],rep(gr,covariateSize))


  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  posCov = orderList$cov[(orderList$cov %in% Reduce(union,H1.all))]
  covariateSizeCar = paste0(covariateSize," covariates")

  # Data to use
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method==buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber==covariateSizeCar,]

  resultCovariateCov <- resultCovariate[resultCovariate$Method==buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]

  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method==buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]

  resultModelCov <- resultModel[resultModel$Method==buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]

  # Value to display
  valueDisplay = list(cov=c(sapply(unique(posCov), function(x){resultCovariateCov[resultCovariateCov$Covariate== x &resultCovariateCov$TypeOfSim=="cov","SelectedInDistinctModel"]*100}),
                            FP = round(mean(resultCovariateCov[!(resultCovariateCov$Covariate %in% unique(posCov))
                                                               & resultCovariateCov$TypeOfSim=="cov",
                                                               "SelectedInDistinctModel"]),digits = 3)*100),

                      corcov=c(sapply(unique(posCov), function(x){resultCovariateCov[resultCovariateCov$Covariate== x
                                                                                     & resultCovariateCov$TypeOfSim=="corcov",
                                                                                     "SelectedInDistinctModel"]*100}),

                               FP = round(mean(resultCovariateCov[!(resultCovariateCov$Covariate %in% unique(posCov))
                                                                  & resultCovariateCov$TypeOfSim=="corcov",
                                                                  "SelectedInDistinctModel"]),digits = 3)*100))



  # Graphs cov
  cov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov",],aes(x=factor(Covariate,levels=orderList$cov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=c(rep(0.5,length(unique(posCov))),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"])) - length(unique(posCov)))),
             color=fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"]))],
             fill = fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"]))])+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With uncorrelated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim=="cov","NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim=="cov","TrueModel"]*100,"%"))+
    geom_segment(x=length(posCov)+0.5,y=0,xend=length(posCov)+0.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title="Criterion :"))+
    annotate("text",x=covariateSize,y=valueDisplay$cov[["FP"]]/100+0.1,hjust=1,color="#9E9FA5",fontface = 'italic',size=5,label=paste0(valueDisplay$cov["FP"],"%"))+
    geom_segment(x=length(posCov)+0.5,y=valueDisplay$cov[["FP"]]/100,xend=covariateSize+0.5,yend=valueDisplay$cov[["FP"]]/100,color="#9E9FA5",linewidth=1,linetype="dashed")+
    coord_cartesian(ylim=c(0,1),clip="off")

  for(i in 1:length(posCov)){
    eval(parse(text='cov <- cov + annotate("text",x=4,y=0.9-(i-1)*0.1,hjust=0,color=fill.vec[i],size=4,
             label=paste0(posCov[i]," :",valueDisplay$cov[[posCov[i]]],"%"))'))
  }

  # Graphs corcov
  corcov =  ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov",],aes(x=factor(Covariate,levels=orderList$corcov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=c(rep(0.5,length(unique(posCov))),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"])) - length(unique(posCov)))),
             color=fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"]))],
             fill = fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"]))])+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With correlated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"]*100,"%"))+
    geom_segment(x=length(posCov)+0.5,y=0,xend=length(posCov)+0.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    theme(legend.position="bottom")+
    guides(fill=guide_legend(title="Criterion :"))+
    annotate("text",x=covariateSize,y=valueDisplay$corcov[["FP"]]/100+0.1,hjust=1,color="#9E9FA5",fontface = 'italic',size=5,label=paste0(valueDisplay$corcov["FP"],"%"))+
    geom_segment(x=length(posCov)+0.5,y=valueDisplay$corcov[["FP"]]/100,xend=covariateSize+0.5,yend=valueDisplay$corcov[["FP"]]/100,color="#9E9FA5",linewidth=1,linetype="dashed")+
    coord_cartesian(ylim=c(0,1),clip="off")

  for(i in 1:length(posCov)){
    eval(parse(text='corcov <- corcov + annotate("text",x=4,y=0.9-(i-1)*0.1,hjust=0,color=fill.vec[i],size=4,
             label=paste0(posCov[i]," :",valueDisplay$corcov[[posCov[i]]],"%"))'))
  }

  # Save plot
  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )

  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection.png"),
           height = 2000, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelection.jpeg"),
           height = 2000, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }


  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelectionCov.png"),
           height = 1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelectionCov.jpeg"),
           height = 1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelectionCorcov.png"),
           height = 1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelectionCorcov.jpeg"),
           height = 1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }
}
