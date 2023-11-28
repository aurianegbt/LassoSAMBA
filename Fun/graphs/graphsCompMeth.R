graphsCompMethod <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files",project,"/H1.all.R"))

  # Color
  gr = "#888888"
  fill.vec = c(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(Reduce(union,H1.all))],rep(gr,covariateSize))


  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = paste0(covariateSize," covariates")
  posCov = orderList$cov[(orderList$cov %in% Reduce(union,H1.all))]

  # Data to use
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method %in% buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber == covariateSizeCar,]

  resultCovariateCov <- resultCovariate[resultCovariate$Method %in% buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]

  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method %in% buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]


  resultModelCov <- resultModel[resultModel$Method %in% buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]

  resultModelCov <- cbind(resultModelCov,text=paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))

  # Value to display
  
  lim=data.frame()
  valueDisplay = data.frame()
  for(t in c("cov","corcov")){
    for(meth in buildMethod){
      lim = rbind(lim,data.frame(TypeOfSim=t,
                                 lim=length(unique(CovariateModelSelectionCov[CovariateModelSelectionCov$Method==meth
                                                                              & CovariateModelSelectionCov$TypeOfSim==t,
                                                                              "Covariate"])),
                                 Method=meth))

      value=c(sapply(posCov, function(y){if(length(resultCovariateCov[resultCovariateCov$TypeOfSim==t &resultCovariateCov$Method == meth &resultCovariateCov$Covariate==y,"SelectedInDistinctModel"])==0){0}else{resultCovariateCov[resultCovariateCov$TypeOfSim==t & resultCovariateCov$Method == meth & resultCovariateCov$Covariate==y,"SelectedInDistinctModel"]}*100}),  FP = round(mean(if(length(resultCovariateCov[resultCovariateCov$TypeOfSim==t & resultCovariateCov$Method == meth & (resultCovariateCov$Covariate %in% unique(unlist(H1.all,use.names = F))),"SelectedInDistinctModel"])==0){0}else{resultCovariateCov[resultCovariateCov$TypeOfSim==t & resultCovariateCov$Method == meth & !(resultCovariateCov$Covariate %in% unique(unlist(H1.all,use.names = F))),"SelectedInDistinctModel"]}*100),digits=2))

      valuemax = max(if(length(resultCovariateCov[resultCovariateCov$TypeOfSim==t & resultCovariateCov$Method == meth & !(resultCovariateCov$Covariate %in% unlist(H1.all,use.names = F)),"SelectedInDistinctModel"])==0){0}else{resultCovariateCov[resultCovariateCov$TypeOfSim==t & resultCovariateCov$Method == meth & !(resultCovariateCov$Covariate %in% unlist(H1.all,use.names = F)),"SelectedInDistinctModel"][(lim[lim$TypeOfSim==t & lim$Method==meth,"lim"]-5):lim[lim$TypeOfSim==t & lim$Method==meth,"lim"]]}*100,na.rm = TRUE)


      covAux = sapply(stringr::str_replace(names(value),"\\.",":"),FUN=function(x){gsub(".*:","",x)},USE.NAMES = FALSE)
      distextAux = stringr::str_c(stringr::str_replace(stringr::str_c(covAux," : "),"FP : ",""),stringr::str_c(unlist(value,use.names = F),"%"))
      coorAux = setNames(rep(0,length(covAux)),covAux)
      coorAux[names(coorAux)=="FP"] <- valuemax
      coorAux[names(coorAux)!="FP"] <- sapply(covAux[covAux!="FP"],function(x){which(posCov==x)},USE.NAMES=F)+0.5
      covAux[covAux!="FP"] <-"cov"

      valueDisplay <- rbind(valueDisplay,data.frame(cov= covAux,
                                                    value= unlist(value,use.names = F),
                                                    type=t,
                                                    Method=meth,
                                                    coor = unname(coorAux),
                                                    distext = distextAux))
    }
  }

  valueModel = cbind(resultModelCov,Parameter=c(names(t.param)[length(H1.all)]),text = paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))


  # Graphs cov
  cov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov",],aes(x=factor(Covariate,levels=orderList$cov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=rep(c(rep(0.5,length(unique(unlist(H1.all,use.names = F)))),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"])) - length(unique(unlist(H1.all,use.names = F))))),length(buildMethod)),
             color=rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"]))],length(buildMethod)),
             fill =rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","Covariate"]))],length(buildMethod)))+
    facet_grid(factor(Method,levels=c("reg","lassoSS","lassoSSCl","lassoSSCrit","elasticnetSS","elasticnetSSCrit","ClustOfVar"))~.,labeller = as_labeller(c(lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step.",lassoSSCrit="Lasso with\nmultiple tresholds.",elasticnetSS="Elastic Net",elasticnetSSCrit="Elastic Net with\nmultiple thresholds.",reg="StepAIC")))+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With uncorrelated covariates, ")+
    geom_segment(x=length(Reduce(union,H1.all))+0.5,y=0,xend=length(Reduce(union,H1.all))+0.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    theme(strip.text.x = element_text(size = 8))+
    ylim(c(0,1))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 14))+
    theme(plot.title = element_text(size=20,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=16))+
    geom_text(data=valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=max(lim[lim$TypeOfSim=="cov","lim"])+0.5,y=(coor+10)/100,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    geom_text(data = resultModelCov[resultModelCov$TypeOfSim=="cov",],mapping=aes(label=text),x=covariateSize+0.5,y=1,hjust=1,size=5,vjust=1)+
    coord_cartesian(ylim=c(0,1),clip="off")

  df = valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="cov",]

  for(i in 1:(nrow(df)/length(buildMethod))){
    eval(parse(text=paste0('cov <- cov + geom_text(data=df[seq(i,nrow(df),nrow(df)/length(buildMethod)),],mapping=aes(label=distext,x=length(Reduce(union,H1.all))+0.6),y=0.9-0.10*(i-1),color=fill.vec[i],size=4,hjust=0)')))
  }

  # Graphs corcov
  corcov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov",],aes(x=factor(Covariate,levels=orderList$corcov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=rep(c(rep(0.5,length(unique(unlist(H1.all,use.names = F)))),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"])) - length(unique(unlist(H1.all,use.names = F))))),length(buildMethod)),
             color=rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"]))],length(buildMethod)),
             fill =rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim=="corcov","Covariate"]))],length(buildMethod)))+
    facet_grid(factor(Method,levels=c("reg","lassoSS","lassoSSCl","lassoSSCrit","elasticnetSS","elasticnetSSCrit","ClustOfVar"))~.,labeller = as_labeller(c(lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step.",lassoSSCrit="Lasso with\nmultiple tresholds.",elasticnetSS="Elastic Net",elasticnetSSCrit="Elastic Net with\nmultiple thresholds.",reg="StepAIC")))+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With correlated covariates, ")+
    geom_segment(x=length(Reduce(union,H1.all))+0.5,y=0,xend=length(Reduce(union,H1.all))+0.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    theme(strip.text.x = element_text(size = 8))+
    ylim(c(0,1))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 14))+
    theme(plot.title = element_text(size=20,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=16))+
    geom_text(data=valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=max(lim[lim$TypeOfSim=="corcov","lim"])+0.5,y=(coor+10)/100,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    geom_text(data = resultModelCov[resultModelCov$TypeOfSim=="corcov",],mapping=aes(label=text),x=covariateSize+0.5,y=1,hjust=1,size=5,vjust=1)+
    coord_cartesian(ylim=c(0,1),clip="off")

  df = valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="cov",]

  for(i in 1:(nrow(df)/length(buildMethod))){
    eval(parse(text=paste0('corcov <- corcov + geom_text(data=df[seq(i,nrow(df),nrow(df)/length(buildMethod)),],mapping=aes(label=distext,x=length(Reduce(union,H1.all))+0.6),y=0.9-0.10*(i-1),color=fill.vec[i],size=4,hjust=0)')))
  }

  # Save plot
  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(stringr::str_wrap(subtitle,100),size=16),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = 3000 + 500*length(buildMethod), width = 4200, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 3000 + 500*length(buildMethod), width = 4200, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )

  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = length(buildMethod)/2*1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = length(buildMethod)/2*1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(stringr::str_wrap(subtitle,100)),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = length(buildMethod)/2*1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = length(buildMethod)/2*1500, width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }
}
