graphsParCompMethod <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Color & covariates
  gr = "#888888"
  fill.vec = c(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(unlist(H1.all,use.names = F))],rep(gr,covariateSize))

  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = paste0(covariateSize," covariates")
  posCov = orderList$cov[(orderList$cov %in% Reduce(union,H1.all))]

  # Change name to display
  newbuildMethod <- c()
  for(k in 1:length(buildMethod)){
    if(buildMethod[k]=="reg"){
      newbuildMethod[k] <- "StepAIC"
    }else if(buildMethod[k]=="lassoSSCl"){
      newbuildMethod[k] <- "Lasso with\nclustering step"
    }else if(buildMethod[k]=="lassoSS"){
      newbuildMethod[k] <- "Lasso"
    }else if(buildMethod[k]=="lasso"){
      newbuildMethod[k] <- "Lasso without\nstability selection"
    }else if(buildMethod[k]=="CoVSS"){
      newbuildMethod[k] <- "ClustOfVar"
    }else if(buildMethod[k]=="elasticnet"){
      newbuildMethod[k] <- "Elastic Net without\nstability selection"
    }else if(buildMethod[k]=="elasticnetSS"){
      newbuildMethod[k] <- "Elastic Net"
    }else if(buildMethod[k]=="lassoSSCrit"){
      newbuildMethod[k] <- "Lasso with\nmultiple thresholds"
    }else if(buildMethod[k]=="elasticnetSSCrit"){
      newbuildMethod[k] <- "Elastic Net with\nmultiple thresholds"
    }
  }

  # Data to use (and change name)
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method %in% buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber == covariateSizeCar,]

  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="reg","Method"] <- "StepAIC"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="lassoSS","Method"] <- "Lasso"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="lasso","Method"] <- "Lasso without\nstability selection"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="lassoSSCrit","Method"] <- "Lasso with\nmultiple thresholds" 
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="elasticnetSS","Method"] <- "Elastic Net"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="elasticnet","Method"] <- "Elastic Net without\nstability selection"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="elasticnetSSCrit","Method"] <- "Elastic Net with\nmultiple thresholds"


  resultCovariateCov <- resultCovariate[resultCovariate$Method %in% buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]

  resultCovariateCov[resultCovariateCov$Method=="reg","Method"] <- "StepAIC"
  resultCovariateCov[resultCovariateCov$Method=="lassoSS","Method"] <- "Lasso"
  resultCovariateCov[resultCovariateCov$Method=="lasso","Method"] <- "Lasso without\nstability selection"
  resultCovariateCov[resultCovariateCov$Method=="lassoSSCrit","Method"] <- "Lasso with\nmultiple thresholds" 
  resultCovariateCov[resultCovariateCov$Method=="elasticnetSS","Method"] <- "Elastic Net"
  resultCovariateCov[resultCovariateCov$Method=="elasticnet","Method"] <- "Elastic Net without\nstability selection"
  resultCovariateCov[resultCovariateCov$Method=="elasticnetSSCrit","Method"] <- "Elastic Net with\nmultiple thresholds"

  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method %in% buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]

  resultCovariateParCov[resultCovariateParCov$Method=="reg","Method"] <- "StepAIC"
  resultCovariateParCov[resultCovariateParCov$Method=="lassoSS","Method"] <- "Lasso"
  resultCovariateParCov[resultCovariateParCov$Method=="lasso","Method"] <- "Lasso without\nstability selection"
  resultCovariateParCov[resultCovariateParCov$Method=="lassoSSCrit","Method"] <- "Lasso with\nmultiple thresholds" 
  resultCovariateParCov[resultCovariateParCov$Method=="elasticnetSS","Method"] <- "Elastic Net"
  resultCovariateParCov[resultCovariateParCov$Method=="elasticnet","Method"] <- "Elastic Net without\nstability selection"
  resultCovariateParCov[resultCovariateParCov$Method=="elasticnetSSCrit","Method"] <- "Elastic Net with\nmultiple thresholds"
  

  resultModelCov <- resultModel[resultModel$Method %in% buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]

  resultModelCov <- cbind(resultModelCov,text=paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))

  resultModelCov[resultModelCov$Method=="reg","Method"] <- "StepAIC"
  resultModelCov[resultModelCov$Method=="lassoSS","Method"] <- "Lasso"
  resultModelCov[resultModelCov$Method=="lasso","Method"] <- "Lasso without\nstability selection"
  resultModelCov[resultModelCov$Method=="lassoSSCrit","Method"] <- "Lasso with\nmultiple thresholds" 
  resultModelCov[resultModelCov$Method=="elasticnetSS","Method"] <- "Elastic Net"
  resultModelCov[resultModelCov$Method=="elasticnet","Method"] <- "Elastic Net without\nstability selection"
  resultModelCov[resultModelCov$Method=="elasticnetSSCrit","Method"] <- "Elastic Net with\nmultiple thresholds"
  
  resultModelParCov <- resultModelPar[resultModelPar$Method %in% buildMethod
                                      & resultModelPar$ProjectNumber==covariateSizeCar,]
  
  resultModelParCov <- cbind(resultModelParCov,text=paste0("Model without False Negatives : ",resultModelParCov$NoFNModel*100,"%"))
  
  resultModelParCov[resultModelParCov$Method=="reg","Method"] <- "StepAIC"
  resultModelParCov[resultModelParCov$Method=="lassoSS","Method"] <- "Lasso"
  resultModelParCov[resultModelParCov$Method=="lasso","Method"] <- "Lasso without\nstability selection"
  resultModelParCov[resultModelParCov$Method=="lassoSSCrit","Method"] <- "Lasso with\nmultiple thresholds" 
  resultModelParCov[resultModelParCov$Method=="elasticnetSS","Method"] <- "Elastic Net"
  resultModelParCov[resultModelParCov$Method=="elasticnet","Method"] <- "Elastic Net without\nstability selection"
  resultModelParCov[resultModelParCov$Method=="elasticnetSSCrit","Method"] <- "Elastic Net with\nmultiple thresholds"
  

  buildMethod <- newbuildMethod

  # Covariate selected
  covariate=list()
  for(t in c("cov","corcov")){
    aux = list()
    for(meth in buildMethod){
      covmetht = setNames(lapply(names(t.param),function(x){unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==t
                                                                                     & CovariateModelSelectionCov$Method==meth
                                                                                     & CovariateModelSelectionCov$Parameter==t.param[x]
                                                                                     ,"Covariate"])}),names(t.param))
      aux <- append(aux,list(covmetht))
      names(aux)[length(aux)] <- meth
    }
    covariate <- append(covariate,list(aux))
    names(covariate)[length(covariate)] <- t
  }

  # To take into account missing covariates for color


  repPar = setNames(lapply(names(t.param),
                     function(p){lapply(setNames(lapply(names(t.param),
                                                  function(p){list(cov= setNames(lapply(buildMethod,
                                                                                  function(m){as.numeric(posCov %in% covariate$cov[[m]][[p]])}),buildMethod),
                                                                   corcov = setNames(lapply(buildMethod,
                                                                                      function(m){as.numeric(posCov %in% covariate$corcov[[m]][[p]])}),buildMethod))}),names(t.param))[[p]],
                                    function(y){lapply(y,
                                                  function(c){c[which(posCov %in% H1.all[[p]])] <- 2;return(c)})})}),names(t.param))

  # Value to Display
  lim = list(cov=length(unique(unlist(covariate$cov))),
             corcov = length(unique(unlist(covariate$cov))))

  valueDisplay = data.frame()
  for(t in c("cov","corcov")){
    for(meth in buildMethod){
      value=setNames(lapply(names(t.param),
                       function(x){c(sapply(H1.all[[x]],
                                      function(y){resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                        & resultCovariateParCov$Method==meth
                                                                        & resultCovariateParCov$TypeOfSim==t
                                                                        &resultCovariateParCov$Covariate==y,"ProportionSelected"]*100}),
                                 FP=round(mean(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                     & resultCovariateParCov$TypeOfSim==t
                                                                     & resultCovariateParCov$Method==meth
                                                                     & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                                     "ProportionSelected"]*100),digits=2))}),names(t.param))

      valuemax = sapply(names(t.param),
                  function(x){max(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                        & resultCovariateParCov$Method==meth
                                                        & resultCovariateParCov$TypeOfSim==t
                                                        & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                        "ProportionSelected"][(lim[[t]]-5):lim[[t]]]*100,na.rm = TRUE)})


      covAux = sapply(stringr::str_replace(names(unlist(value)),"\\.",":"),FUN=function(x){gsub(".*:","",x)},USE.NAMES = FALSE)
      distextAux = stringr::str_c(stringr::str_replace(stringr::str_c(covAux," : "),"FP : ",""),
                                  stringr::str_c(unlist(value,use.names = F),"%"))
      coorAux = setNames(rep(0,length(covAux)),covAux)
      coorAux[names(coorAux)=="FP"] <- valuemax
      coorAux[names(coorAux)!="FP"] <- sapply(covAux[covAux!="FP"],function(x){which(posCov==x)},USE.NAMES=F)+0.6
      covAux[covAux!="FP"] <-"cov"


      repAuxPar = setNames(rep(1,length(names(t.param))),names(t.param))
      repAuxPar[names(H1.all)] <- sapply(H1.all,length)+1

      valueDisplay <- rbind(valueDisplay,data.frame(cov= covAux,
                                                    value= unlist(value,use.names = F),
                                                    type=t,
                                                    Method=meth,
                                                    coor = unname(coorAux),
                                                    Parameter= unname(rep(t.param,unname(repAuxPar))),
                                                    distext = distextAux))
    }
  }
  valueModel = cbind(resultModelParCov,Parameter=rev(t.param)[1])


  # Cov graphs
  precious = rep(c("#468b97","#ef6262", "#74C385","#8e6aa0","#ee6c4d","#007194")[1:length(unlist(H1.all,use.names = F))],length(buildMethod))

  cmdA = list()
  cmdC = list()
  covPast = 0
  colDisplay=c()
  for(meth in buildMethod){
    for(p in rev(names(t.param))){
      cAux <- repPar[[p]]$cov[[meth]]
      if(length(cAux[cAux!=0])==0){
        cmdA <- append(cmdA,paste0('rep(0.25,length(covariate$cov[["',meth,'"]][["',p,'"]]))'))
      }else{
        cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariate$cov[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))

      }

      if(length(cAux[cAux!=0])==0){
        cmdC <- append(cmdC, paste0('rep(gr,length(covariate$cov[["',meth,'"]][["',p,'"]]))'))
      }else{
        if(length(cAux[cAux==2])!=0){
          colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
        }
        cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
        cAux[cAux==1] <- gr
        cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariate$cov[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))
      }

      covPast = covPast + sum(cAux %in% precious)
    }
  }
  eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
  eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))


  cov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim=="cov",],aes(x=factor(Covariate,levels=orderList$cov)))+
    geom_bar(position=position_dodge(preserve = "single"),alpha=alphaScale,
             color= colorScale,
             fill = colorScale)+
    facet_nested(factor(Method,levels=c("StepAIC","Lasso","Elastic Net","Lasso with\nclustering step","Lasso with\nmultiple thresholds","Elastic Net with\nclustering step","Elastic Net with\nmultiple thresholds"))+Parameter~.,labeller=labeller(Parameter=label_parsed))+
    xlab("Covariates")+
    ylab("Count")+
    ggtitle("With uncorrelated covariates, ")+
    theme(axis.text.x = element_text(size = 10, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    ylim(c(0,max(resultCovariateCov[resultCovariateCov$TypeOfSim=="cov","NumberofModel"])))+
    theme(axis.title = element_text(size=20))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=25,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])*(0.90-(coor-2)*0.2),hjust=0),color=rev(colDisplay),size=7)+
    geom_text(data=valueDisplay[valueDisplay$type=="cov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$cov+0.5,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=7)+
    geom_text(data=valueModel[valueModel$TypeOfSim=="cov",],x=lim$cov+0.5,y=max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"]),hjust=1,size=7,vjust=1,mapping=aes(label=text))+
    coord_cartesian(ylim=c(0,max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])),clip="off")


  # Corcov graphs
  cmdA = list()
  cmdC = list()
  covPast = 0
  colDisplay=c()
  for(meth in buildMethod){
    for(p in rev(names(t.param))){
      cAux <- repPar[[p]]$corcov[[meth]]
      if(length(cAux[cAux!=0])==0){
        cmdA <- append(cmdA,paste0('rep(0.25,length(covariate$corcov[["',meth,'"]][["',p,'"]]))'))
      }else{
        cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariate$corcov[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))

      }

      if(length(cAux[cAux!=0])==0){
        cmdC <- append(cmdC, paste0('rep(gr,length(covariate$corcov[["',meth,'"]][["',p,'"]]))'))
      }else{
        if(length(cAux[cAux==2])!=0){
          colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
        }
        cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
        cAux[cAux==1] <- gr
        cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariate$corcov[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))
      }

      covPast = covPast + sum(cAux %in% precious)
    }
  }
  eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
  eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))

  corcov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim=="corcov",],aes(x=factor(Covariate,levels=orderList$corcov)))+
    geom_bar(position=position_dodge(preserve = "single"),alpha=alphaScale,
             color= colorScale,
             fill = colorScale)+
    facet_nested(factor(Method,levels=c("StepAIC","Lasso","Elastic Net","Lasso with\nclustering step","Lasso with\nmultiple thresholds","Elastic Net with\nclustering step","Elastic Net with\nmultiple thresholds"))+Parameter~.,labeller=labeller(Parameter=label_parsed))+
    xlab("Covariates")+
    ylab("Count")+
    ggtitle("With correlated covariates, ")+
    theme(axis.text.x = element_text(size = 10, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    ylim(c(0,max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])))+
    theme(axis.title = element_text(size=18))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=25,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])*(0.90-(coor-2)*0.2),hjust=0),color=rev(colDisplay),size=7)+
    geom_text(data=valueDisplay[valueDisplay$type=="corcov" & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$corcov+0.5,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=7)+
    geom_text(data=valueModel[valueModel$TypeOfSim=="corcov",],x=lim$corcov+0.5,y=max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"]),hjust=1,size=7,vjust=1,mapping=aes(label=text))+
    coord_cartesian(ylim=c(0,max(resultCovariateParCov[resultCovariateParCov$TypeOfSim=="cov","NumberofModel"])),clip="off") 
  
  
  # Save plot
  annotate_figure(ggarrange(cov, corcov, nrow = 2),
                  top=text_grob(stringr::str_wrap(subtitle,80),size=25),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=20,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=30,color="#862B0D")
    )



  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = length(buildMethod)*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE, bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
            height = length(buildMethod)*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE,device=grDevices::jpeg)
  }

  annotate_figure(cov,
                  top=text_grob(stringr::str_wrap(subtitle,80),size=25),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=20,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=30,color="#862B0D")
    )


  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
            height = length(buildMethod)/2*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE, bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
            height = length(buildMethod)/2*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE,device=grDevices::jpeg)
  }

  annotate_figure(corcov,
                  top=text_grob(stringr::str_wrap(subtitle,80),size=25),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=20,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=30,color="#862B0D")
    )


  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
            height = length(buildMethod)/2*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE, bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
            height = length(buildMethod)/2*list("3"=2500,"4"=3000,"5"=3500)[[as.character(length(t.param))]], width = list("10"=4500,"50"=5000,"200"=6000,"500"=8000)[[as.character(covariateSize)]], units="px",limitsize =   FALSE,device=grDevices::jpeg)
  }
}
