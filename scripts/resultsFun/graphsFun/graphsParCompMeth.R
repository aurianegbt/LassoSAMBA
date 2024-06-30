graphsParCompMethod <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Color & covariates
  gr = "#888888"
  fill.vec = c(c("#468b97", "#ef6262", "#74C385", "#8e6aa0", "#ee6c4d",                   "#1e90ff", "#ffa500", "#ff69b4", "#32cd32", "#4169e1",                  "#ff6347", "#6a5acd", "#20b2aa", "#f08080", "#6495ed",                  "#9acd32", "#9370db", "#00ced1", "#ff4500", "#7b68ee",                  "#2e8b57", "#ba55d3", "#00bfff", "#d2691e", "#4682b4")[1:length(unlist(H1.all,use.names = F))],rep(gr,200))

  posCov = orderList[(orderList %in% Reduce(union,H1.all))]

  # Change name to display
  newbuildMethod <- c()
  for(k in 1:length(buildMethod)){
    if(buildMethod[k]=="reg"){
      newbuildMethod[k] <- "StepAIC\nwith stat. test"
    }else if(buildMethod[k]=="lasso"){
      newbuildMethod[k] <- "Lasso without\nstability selection"
    }else if(buildMethod[k]=="elasticnet"){
      newbuildMethod[k] <- "Elastic Net without\nstability selection"
    }else if(buildMethod[k]=="lassoSS"){
      newbuildMethod[k] <- "Lasso\nwith stat. test"
    }else if(buildMethod[k]=="elasticnetSS"){
      newbuildMethod[k] <- "Elastic Net\nwith stat. test"
    }else if(stringr::str_detect(buildMethod[k],"regPEN")){
      newbuildMethod[k] <- paste0("penalized stepAIC\npen=",stringr::str_remove(buildMethod[k],"regPEN"))
    }else if(buildMethod[k]=="lassoSSnoCov0"){
      newbuildMethod[k] <- "Lasso"
    }else if(buildMethod[k]=="lassoSSCritnoCov0"){
      newbuildMethod[k] <- "Lasso\nwith adapted thresholds"
    }else if(buildMethod[k]=="regnoCov0"){
      newbuildMethod[k] <- "StepAIC"
    }else if(buildMethod[k]=="elasticnetSSnoCov0"){
      newbuildMethod[k] <- "Elastic Net"
    }else if(buildMethod[k]=="sharp"){
      newbuildMethod[k] <- "Lasso calibrated using sharp\nwith stat. test"
    }else if(buildMethod[k]=="sharpnoCov0"){
      newbuildMethod[k] <- "Lasso calibrated using sharp"
    }else if(stringr::str_detect(buildMethod[k],"sharpnoCov0FDP")){
      newbuildMethod[k] <- paste0("Lasso using stability score\nE[FDR]<",stringr::str_remove(buildMethod[k],"sharpnoCov0FDP"),"%")
    }else if(stringr::str_detect(buildMethod[k],"sharpnoCov0")){
      newbuildMethod[k] <- paste0("Lasso using stability score\n",stringr::str_remove(buildMethod[k],"sharpnoCov0"),"% higher score")
    }else if(stringr::str_detect(buildMethod[k],"sharpnoCov0FDP")){
      newbuildMethod[k] <- paste0("Lasso using stability score\nE[FDR]<",stringr::str_remove(buildMethod[k],"sharpnoCov0FDP"),"%")
    }else if(stringr::str_detect(buildMethod[k],"sharp")){
      newbuildMethod[k] <- paste0("Lasso using stability score\n",stringr::str_remove(buildMethod[k],"sharp"),"% higher score\nwith stat. test")
    }
  }
  
  # Data to use (and change name)
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method %in% buildMethod,]

  resultCovariateCov <- resultCovariate[resultCovariate$Method %in% buildMethod,]
  
  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method %in% buildMethod,]
  
  resultModelCov <- resultModel[resultModel$Method %in% buildMethod,]

  resultModelCov <- cbind(resultModelCov,text=paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))
  
  resultModelParCov <- resultModelPar[resultModelPar$Method %in% buildMethod,]
  
  resultModelParCov <- cbind(resultModelParCov,text=paste0("Model without False Negatives : ",resultModelParCov$NoFNModel*100,"%"))
  
  for(meth in buildMethod){       
    CovariateModelSelectionCov[CovariateModelSelectionCov$Method==meth,"Method"] <- newbuildMethod[buildMethod==meth]
   resultCovariateCov[resultCovariateCov$Method==meth,"Method"] <- newbuildMethod[buildMethod==meth]
   resultCovariateParCov[resultCovariateParCov$Method==meth,"Method"] <- newbuildMethod[buildMethod==meth]
   resultModelCov[resultModelCov$Method==meth,"Method"] <- newbuildMethod[buildMethod==meth]
   resultModelParCov[resultModelParCov$Method==meth,"Method"] <- newbuildMethod[buildMethod==meth]
  }

  buildMethod <- newbuildMethod

  # Covariate selected
  covariate = list()
  for(meth in buildMethod){
    covmetht = setNames(lapply(names(t.param),function(x){unique(CovariateModelSelectionCov[ CovariateModelSelectionCov$Method==meth & CovariateModelSelectionCov$Parameter==t.param[x],"Covariate"])}),
                        names(t.param))
    covariate <- append(covariate,list(covmetht))
    names(covariate)[length(covariate)] <- meth
  }
  
  # To take into account missing covariates for color
  repPar = setNames(lapply(names(t.param),
                           function(p){setNames(lapply(buildMethod,FUN=function(m){as.numeric(posCov %in% covariate[[m]][[p]])}),buildMethod)})
                    ,names(t.param))
  repPar = setNames(lapply(names(t.param),
                           function(p){setNames(lapply(repPar[[p]],FUN=function(rep){rep[which(posCov %in% H1.all[[p]])] <-2;return(rep)}),buildMethod)}),names(t.param))

  # Value to Display
  lim = length(unique(unlist(covariate)))

  valueDisplay = data.frame()
  for(meth in buildMethod){
    value=setNames(lapply(names(t.param),
                          function(x){c(sapply(H1.all[[x]],
                                               function(y){resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                                 & resultCovariateParCov$Method==meth
                                                                                 & resultCovariateParCov$Covariate==y,"ProportionSelected"]*100}),
                                        meanSelFP=round(mean(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                                   & resultCovariateParCov$Method==meth
                                                                                   & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                                                   "ProportionSelected"]*100),digits=2))}),names(t.param))
    
    valuemax = sapply(names(t.param),
                      function(x){max(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                            & resultCovariateParCov$Method==meth
                                                            & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                            "ProportionSelected"][(lim-5):lim]*100,na.rm = TRUE)})
    
    
    covAux = sapply(stringr::str_replace(names(unlist(value)),"\\.",":"),FUN=function(x){gsub(".*:","",x)},USE.NAMES = FALSE)
    distextAux = stringr::str_c(stringr::str_replace(stringr::str_c(covAux," : "),"meanSelFP : ",""),
                                stringr::str_c(unlist(value,use.names = F),"%"))
    coorAux = setNames(rep(0,length(covAux)),covAux)
    coorAux[names(coorAux)=="meanSelFP"] <- valuemax
    coorAux[names(coorAux)!="meanSelFP"] <- sapply(covAux[covAux!="meanSelFP"],function(x){which(posCov==x)},USE.NAMES=F)+0.6
    covAux[covAux!="meanSelFP"] <-"cov"
    
    
    repAuxPar = setNames(rep(1,length(names(t.param))),names(t.param))
    repAuxPar[names(H1.all)] <- sapply(H1.all,length)+1
    
    valueDisplay <- rbind(valueDisplay,data.frame(cov= covAux,
                                                  value= unlist(value,use.names = F),
                                                  Method=meth,
                                                  coor = unname(coorAux),
                                                  Parameter= unname(rep(t.param,unname(repAuxPar))),
                                                  distext = distextAux))
  }
  valueDisplay$coor2 <- rep(max(valueDisplay[valueDisplay$cov=="cov",'coor']),nrow(valueDisplay[valueDisplay$cov=="cov",] ))
  valueModel = cbind(resultModelParCov,Parameter=rev(t.param)[1])


  precious = rev(rep(c("#468b97", "#ef6262", "#74C385", "#8e6aa0", "#ee6c4d",                   "#1e90ff", "#ffa500", "#ff69b4", "#32cd32", "#4169e1",                  "#ff6347", "#6a5acd", "#20b2aa", "#f08080", "#6495ed",                  "#9acd32", "#9370db", "#00ced1", "#ff4500", "#7b68ee",                  "#2e8b57", "#ba55d3", "#00bfff", "#d2691e", "#4682b4")[1:length(unlist(H1.all,use.names = F))],length(buildMethod)))

  cmdA = list()
  cmdC = list()
  covPast = 0
  colDisplay=c()
  for(meth in buildMethod){
    for(p in rev(names(t.param))){
      cAux <- repPar[[p]][[meth]]
      if(length(cAux[cAux!=0])==0){
        cmdA <- append(cmdA,paste0('rep(0.25,length(covariate[["',meth,'"]][["',p,'"]]))'))
      }else{
        cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariate[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))

      }

      if(length(cAux[cAux!=0])==0){
        cmdC <- append(cmdC, paste0('rep(gr,length(covariate[["',meth,'"]][["',p,'"]]))'))
      }else{
        if(length(cAux[cAux==2])!=0){
          colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
        }
        cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
        cAux[cAux==1] <- gr
        cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariate[["',meth,'"]][["',p,'"]]) - ',sum(cAux!=0),'))'))
      }

      covPast = covPast + sum(cAux %in% precious)
    }
  }
  eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
  eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))

  plot = ggplot(CovariateModelSelectionCov,aes(x=factor(Covariate,levels=orderList)))+
    geom_bar(position=position_dodge(preserve = "single"),alpha=alphaScale,
             color= colorScale,
             fill = colorScale)+
    facet_nested(factor(Method,levels=newbuildMethod)+Parameter~.,labeller=labeller(Parameter=label_parsed))+
    xlab("Covariates")+
    ylab("Count")+
    theme(axis.text.x = element_text(size = 8, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    ylim(c(0,lim))+
    theme(axis.title = element_text(size=20))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=25,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[ valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor2,y=if(project=="Pasin" | project=="GaussianPasin"){max(resultCovariateParCov[,"NumberofModel"])*0.9}else{max(resultCovariateParCov[,"NumberofModel"])*(0.80-(coor-2)*0.4)},hjust=0,vjust=1),color=rev(colDisplay),size=7)+
    geom_text(data=valueDisplay[ valueDisplay=="meanSelFP",],mapping=aes(label=distext,x=+Inf,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=7)+
    coord_cartesian(ylim=c(0,max(resultCovariateParCov[,"NumberofModel"])),clip="off")
  
  
  # Save plot
  annotate_figure(plot,
                  top=text_grob(stringr::str_wrap(subtitle,120),size=18),
  ) %>%
    annotate_figure(
      top=text_grob("covariates presence in final model with parameters link",
                    face="italic",size=20,color="#9c9c9c")
    )%>%
    annotate_figure(
      top=text_grob("Covariate Selection Frequency",
                    face="bold",size=25,color="#862B0D")
    )



  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter.png"),
           height = length(buildMethod)*1000, width = 5000, units="px",limitsize =   FALSE, bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter.jpeg"),
            height = length(buildMethod)*1000, width = 5000, units="px",limitsize =   FALSE,device=grDevices::jpeg)
  }
}
