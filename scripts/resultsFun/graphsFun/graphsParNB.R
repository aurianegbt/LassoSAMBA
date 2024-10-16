graphsParNB <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){
  
  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))
  
  # Color & covariates
  gr = "#888888"
  fill.vec = c(c("#468b97", "#ef6262", "#74C385", "#8e6aa0", "#ee6c4d", 
                 "#1e90ff", "#ffa500", "#ff69b4", "#32cd32", "#4169e1",
                 "#ff6347", "#6a5acd", "#20b2aa", "#f08080", "#6495ed",
                 "#9acd32", "#9370db", "#00ced1", "#ff4500", "#7b68ee",
                 "#2e8b57", "#ba55d3", "#00bfff", "#d2691e", "#4682b4")
[1:length(unlist(H1.all,use.names = F))],rep(gr,200))
  
  posCov = orderList[(orderList %in% Reduce(union,H1.all))]
  
  # Data.frame to use
  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method==buildMethod,]
  
  if(nrow(CovariateModelSelectionCov)==0){
    return(NULL)
  }else{
    
    resultCovariateCov <- resultCovariate[resultCovariate$Method==buildMethod,]
    
    resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method==buildMethod,]
    
    resultModelCov <- resultModel[resultModel$Method==buildMethod,]
    
    errorStatsParCov <- errorStatsPar[errorStatsPar$Method==buildMethod,]
    
    covariatePar = setNames(lapply(names(t.param),
                                   function(x){
                                     unique(CovariateModelSelectionCov[CovariateModelSelectionCov$Parameter==t.param[x],"Covariate"])}),names(t.param))
    
    
    lim = length(Reduce(union,covariatePar))
    
    
    # To take into account missing covariates for color
    repPar = setNames(lapply(names(t.param),
                             function(p){as.numeric(posCov %in% covariatePar[[p]])})
                      ,names(t.param))
    repPar = setNames(lapply(names(t.param),
                             function(p){repPar[[p]][which(posCov %in% H1.all[[p]])] <-2;return(repPar[[p]])}),names(t.param))
    
    
    
    ## Value to display
    value = setNames(lapply(names(t.param),
                            function(x){c(sapply(H1.all[[x]],
                                                 function(y){resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                                   & resultCovariateParCov$Covariate==y,"ProportionSelected"]*100}),
                                          meanSelFP=round(mean(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                                                     & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                                                     "ProportionSelected"]*100),digits=2))}),names(t.param))
    
    
    valuemax = sapply(names(t.param),
                      function(x){max(resultCovariateParCov[resultCovariateParCov$Parameter==t.param[x]
                                                            & !(resultCovariateParCov$Covariate %in% H1.all[[x]]),
                                                            "ProportionSelected"]*100,na.rm = TRUE)})
    
    
    covAux = sapply(stringr::str_replace(names(unlist(value)),"\\.",":"),FUN=function(x){gsub(".*:","",x)},USE.NAMES = FALSE)
    distextAux = stringr::str_c(stringr::str_replace(stringr::str_c(covAux," : "),"meanSelFP : ",""),
                                stringr::str_c(unlist(value,use.names = F),"%"))
    coorAux = setNames(rep(0,length(covAux)),covAux)
    coorAux[names(coorAux)=="meanSelFP"] <- valuemax
    coorAux[names(coorAux)!="meanSelFP"] <- sapply(names(value),function(p){
      rev(1:length(H1.all[[p]]))+0.5
    },USE.NAMES=F)
    covAux[covAux!="meanSelFP"] <- "cov"
    
    coor2Aux <- coorAux
    coor2Aux[names(coorAux)!="meanSelFP"] <- sapply(names(value),FUN=function(p){
      rep(which(posCov==H1.all[[p]][length(H1.all[[p]])]),length(H1.all[[p]]))+0.5
    })
    
    repAuxPar = setNames(rep(1,length(names(t.param))),names(t.param))
    repAuxPar[names(H1.all)] <- sapply(H1.all,length)+1
    
    valueDisplay <- data.frame(cov= covAux,
                               value= unlist(value,use.names = F),
                               coor = unname(coorAux),
                               coor2 = unname(coor2Aux),
                               Parameter= unname(rep(t.param,unname(repAuxPar))),
                               distext = distextAux)
    
    # Color and fill arguments / cov Graphs
    precious = rev(c("#468b97", "#ef6262", "#74C385", "#8e6aa0", "#ee6c4d",                   "#1e90ff", "#ffa500", "#ff69b4", "#32cd32", "#4169e1",                  "#ff6347", "#6a5acd", "#20b2aa", "#f08080", "#6495ed",                  "#9acd32", "#9370db", "#00ced1", "#ff4500", "#7b68ee",                  "#2e8b57", "#ba55d3", "#00bfff", "#d2691e", "#4682b4")[1:length(unlist(H1.all,use.names = F))])
    cmdA = list()
    cmdC = list()
    covPast = 0
    colDisplay=c()
    for(p in rev(names(t.param))){
      cAux <- repPar[[p]]
      if(length(cAux[cAux!=0])==0){
        cmdA <- append(cmdA,paste0('rep(0.25,length(covariatePar[["',p,'"]]))'))
      }else{
        cmdA <- append(cmdA,paste0('c(',paste0(0.25*cAux[cAux!=0],collapse=","),',rep(0.25,length(covariatePar[["',p,'"]]) - ',sum(cAux!=0),'))'))
        
      }
      
      if(length(cAux[cAux!=0])==0){
        cmdC <- append(cmdC, paste0('rep(gr,length(covariatePar[["',p,'"]]))'))
      }else{
        if(length(cAux[cAux==2])!=0){
          colDisplay <- c(colDisplay,rev(precious[(1+covPast):(covPast+length(cAux[cAux==2]))]))
        }
        cAux[cAux==2] <- precious[(1+covPast):(covPast+length(cAux[cAux==2]))]
        cAux[cAux==1] <- gr
        cmdC <- append(cmdC, paste0('c("',paste0(cAux[cAux!=0],collapse='","'),'",rep(gr,length(covariatePar[["',p,'"]]) - ',sum(cAux!=0),'))'))
      }
      
      covPast = covPast + sum(cAux %in% precious)
    }
    eval(parse(text=paste0('alphaScale = c(',paste0(cmdA,collapse=","),')')))
    eval(parse(text=paste0('colorScale = c(',paste0(cmdC,collapse=","),')')))
    
    plot =
      ggplot(CovariateModelSelectionCov,aes(x=factor(Covariate,levels=orderList)))+
      geom_bar(position=position_dodge(preserve = "single"), alpha=alphaScale,color=colorScale, fill =colorScale)+
      facet_grid(Parameter~.,labeller=label_parsed)+
      xlab("Covariates")+
      ylab("Count")+
      ggtitle(" ",
              subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[,"NoFNModel"]*100,"%","\n",
                               "   • Final model is the true one  : ",resultModelCov[,"TrueModel"]*100,"%"))+
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
      geom_text(data=valueDisplay[ valueDisplay$cov=="meanSelFP",],mapping=aes(label=distext,x=+Inf,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
      coord_cartesian(ylim=c(0,unique(resultCovariateParCov[,"NumberofModel"])),clip="off")+
      geom_text(data=valueDisplay[ valueDisplay$cov=="cov",],mapping=aes(x=coor2,label=distext,y=if(project=="Pasin"){unique(resultCovariateParCov[,"NumberofModel"])}else{rev(unique(resultCovariateParCov[,"NumberofModel"])*(0.80-(coor-2)*0.3))},hjust=0,vjust=1),color=rev(colDisplay),size=5)
    
    # Save plot
    gp <- 
      annotate_figure(plot,
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
      ggsave(plot = gp,filename = paste0(Folder,"/NumberSelectionParameter.png"),
             height = 1500,width = 2500, units = "px", bg='transparent',device=grDevices::png)
    }
    if(JPEG){
      ggsave(plot=gp, filename = paste0(Folder,"/NumberSelectionParameter.jpeg"),
             height = 1500,width = 2500, units = "px",device=grDevices::jpeg)
    }
  }
}