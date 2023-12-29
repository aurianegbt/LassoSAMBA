graphsParNB <- function(Folder,subtitle,project,covariateSize, buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Color & covariates
  gr = "#888888"
  fill.vec = c(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(unlist(H1.all,use.names = F))],rep(gr,covariateSize))

  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  posCov = orderList$cov[(orderList$cov %in% Reduce(union,H1.all))]
  covariateSizeCar = paste0(covariateSize," covariates")

  # Data.frame to use
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method==buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber == covariateSizeCar,]

  resultCovariateCov <- resultCovariate[resultCovariate$Method==buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]

  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method==buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]

  resultModelCov <- resultModel[resultModel$Method==buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]

  covariatePar = setNames(lapply(names(t.param),
                        function(x){list(
                          cov = unique(CovariateModelSelectionCov[
                            CovariateModelSelectionCov$TypeOfSim=="cov"
                            & CovariateModelSelectionCov$Parameter==t.param[x],"Covariate"]),
                          corcov = unique(CovariateModelSelectionCov[
                            CovariateModelSelectionCov$TypeOfSim=="corcov"
                            & CovariateModelSelectionCov$Parameter==t.param[x],"Covariate"]))}),names(t.param))


  lim = list(cov=length(Reduce(union,lapply(covariatePar,function(x){x$cov}))),
             corcov = length(Reduce(union,lapply(covariatePar,function(x){x$corcov}))))


  # To take into account missing covariates for color
  repPar = setNames(lapply(names(t.param),
                           function(x){lapply(setNames(lapply(names(t.param),
                                                         function(p){list(
                                                          cov= as.numeric(posCov %in% covariatePar[[p]]$cov),
                                                          corcov = as.numeric(posCov %in% covariatePar[[p]]$corcov))}),
                                                  names(t.param))[[x]],
                                           function(y){y[which(posCov %in% H1.all[[x]])] <- 2;return(y)})}),
                    names(t.param))



  ## Value to display
  valueDisplay = data.frame()
  for(t in c("cov","corcov")){
    value = setNames(lapply(names(t.param),
                      function(x){c(sapply(H1.all[[x]],
                                     function(y){resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                       & resultCovariateParCov$TypeOfSim==t
                                                                       & resultCovariateParCov$Covariate==y,"ProportionSelected"]*100}),
                                    FP=round(mean(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                        & resultCovariateParCov$TypeOfSim==t
                                                                        & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                                        "ProportionSelected"]*100),digits=2))}),names(t.param))


    valuemax = sapply(names(t.param),
                function(x){max(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                      & resultCovariateParCov$TypeOfSim==t
                                                      & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                      "ProportionSelected"][(lim[[t]]-5):lim[[t]]]*100,na.rm = TRUE)})


    covAux = sapply(stringr::str_replace(names(unlist(value)),"\\.",":"),FUN=function(x){gsub(".*:","",x)},USE.NAMES = FALSE)
    distextAux = stringr::str_c(stringr::str_replace(stringr::str_c(covAux," : "),"FP : ",""),
                                stringr::str_c(unlist(value,use.names = F),"%"))
    coorAux = setNames(rep(0,length(covAux)),covAux)
    coorAux[names(coorAux)=="FP"] <- valuemax
    coorAux[names(coorAux)!="FP"] <- sapply(covAux[covAux!="FP"],function(x){which(posCov==x)},USE.NAMES=F)+0.5
    covAux[covAux!="FP"] <-"cov"


    repAuxPar = setNames(rep(1,length(names(t.param))),names(t.param))
    repAuxPar[names(H1.all)] <- sapply(H1.all,length)+1

    valueDisplay <- rbind(valueDisplay,data.frame(cov= covAux,
                                                  value= unlist(value,use.names = F),
                                                  type=t,
                                                  coor = unname(coorAux),
                                                  Parameter= unname(rep(t.param,unname(repAuxPar))),
                                                  distext = distextAux))

  }

  # Color and fill arguments / cov Graphs
  precious = rev(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(unlist(H1.all,use.names = F))])
  cmdA = list()
  cmdC = list()
  covPast = 0
  colDisplay=c()
  for(p in rev(names(t.param))){
    cAux <- repPar[[p]]$cov
    if(length(cAux[cAux!=0])==0){
      cmdA <- append(cmdA,paste0('rep(0.25,length(covariatePar[["',p,'"]]$cov))'))
    }else{
      cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariatePar[["',p,'"]]$cov) - ',sum(cAux!=0),'))'))

    }

    if(length(cAux[cAux!=0])==0){
      cmdC <- append(cmdC, paste0('rep(gr,length(covariatePar[["',p,'"]]$cov))'))
    }else{
      if(length(cAux[cAux==2])!=0){
        colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
      }
      cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
      cAux[cAux==1] <- gr
      cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariatePar[["',p,'"]]$cov) - ',sum(cAux!=0),'))'))
    }

    covPast = covPast + sum(cAux %in% precious)
  }
  eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
  eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))


  cov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim=="cov",],aes(x=factor(Covariate,levels=orderList$cov)))+
    geom_bar(position=position_dodge(preserve = "single"), alpha=alphaScale,color=colorScale, fill =colorScale)+
    facet_grid(Parameter~.,labeller=label_parsed)+
    xlab("Covariates")+
    ylab("Count")+
    ggtitle("With uncorrelated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim=="cov","NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim=="cov","TrueModel"]*100,"%"))+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    ylim(c(0,unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$cov,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    coord_cartesian(ylim=c(0,unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])),clip="off")+
    geom_text(data=valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=if(project=="Pasin"){unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])}else{rev(unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])*(0.90-(coor-2)*0.1)-2)},hjust=0,vjust=1),color=rev(colDisplay),size=5)

  # Color and fill arguments / corcov Graphs
  cmdA = list()
  cmdC = list()
  covPast = 0
  colDisplay=c()
  for(p in rev(names(t.param))){
    cAux <- repPar[[p]]$corcov
    if(length(cAux[cAux!=0])==0){
      cmdA <- append(cmdA,paste0('rep(0.25,length(covariatePar[["',p,'"]]$corcov))'))
    }else{
      cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariatePar[["',p,'"]]$corcov) - ',sum(cAux!=0),'))'))

    }

    if(length(cAux[cAux!=0])==0){
      cmdC <- append(cmdC, paste0('rep(gr,length(covariatePar[["',p,'"]]$corcov))'))
    }else{
      if(length(cAux[cAux==2])!=0){
        colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
      }
      cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
      cAux[cAux==1] <- gr
      cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariatePar[["',p,'"]]$corcov) - ',sum(cAux!=0),'))'))
    }

    covPast = covPast + sum(cAux %in% precious)
  }
  eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
  eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))


  corcov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim=="corcov",],aes(x=factor(Covariate,levels=orderList$corcov)))+
    geom_bar(position=position_dodge(preserve = "single"), alpha=alphaScale,color=colorScale, fill =colorScale)+
    facet_grid(Parameter~.,labeller=label_parsed)+
    xlab("Covariates")+
    ylab("Count")+
    ggtitle("With correlated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"]*100,"%"))+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    ylim(c(0,unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="corcov","NumberofModel"])))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$corcov,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    coord_cartesian(ylim=c(0,unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])),clip="off")+
    geom_text(data=valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=if(project=="Pasin"){unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="corcov","NumberofModel"])}else{rev(unique(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="corcov","NumberofModel"])*(0.90-(coor-2)*0.1)-2)},hjust=0,vjust=1),color=rev(colDisplay),size=5)

  # Save plot
  annotate_figure(ggarrange(cov, corcov, nrow = 2),
                  top=text_grob(stringr::str_wrap(subtitle,100)),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=12,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=18,color="#862B0D")
    )

  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter.png"),
           height = 2000+400*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter.jpeg"),
           height = 2000+400*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }

  annotate_figure(cov,
                  top=text_grob(stringr::str_wrap(subtitle,100)),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=12,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=18,color="#862B0D")
    )
  
  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameterCov.png"),
           height = 1000+200*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameterCov.jpeg"),
           height = 1000+200*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }

annotate_figure(corcov,
                top=text_grob(stringr::str_wrap(subtitle,100)),
) %>%
  annotate_figure(
    top=text_grob("covariates presence in final model with parameters link",
                  face="italic",size=12,color="#9c9c9c")
  )%>%
  annotate_figure(
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=18,color="#862B0D")
  )

  
  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameterCorcov.png"),
           height =1000+200*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameterCorcov.jpeg"),
           height =1000+200*length(t.param),width = list("10"=2500,"50"=3000,"200"=4000,"500"=5000)[[as.character(covariateSize)]], units = "px",device=grDevices::jpeg)
  }


}
