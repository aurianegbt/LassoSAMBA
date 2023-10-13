graphsGenerate <- function(covariateSize=200,
                           buildMethod=c("regression","lassoSS"),
                           buildOption="standard",
                           Rsmlx="",
                           JPEG = FALSE,
                           PNG = TRUE){
  if(file.exists(paste0("Results/Results",Rsmlx))){
    compTime <- as.logical(readLines(paste0("Results/Results",Rsmlx,"/compTime.txt")))
  }else{
    compTime = FALSE
  }

  library(ggplot2,quietly=TRUE)
  library(ggpubr,quietly=TRUE)
  library(webshot2)
  suppressMessages(library(flextable))
  library(grid)
  library(gtable)
  library(gridExtra)
  library(scales)
  source("~/Travail/00_Theme.R")
  source("Fun/graphs.R")


  generalsubtitle = paste0("Among 100 simulated datasets of simulated humoral immune response to a Ebola Vaccine (Ad26/MVA) for 100 individuals, with ",covariateSize," covariates.\n")
  Titlelist = list(regression = list(Mstar = "Model built with SAMBA.\nThe Model is initialized with the true one.",
                                     noCov0="Model built with SAMBA.\nThe Model is build whithout excluding covariates.",
                                     standard ="Model built with SAMBA."),
                   lassoSS = list(Mstar = "Model built with a lasso approach within SAMBA.\nThe Model is initialized with the true one.",
                                  noCov0="Model built with a lasso approach within SAMBA.\nThe Model is build whithout excluding covariates.",
                                  standard ="Model built with a lasso approach within SAMBA."),
                   lassoSSCl = list(Mstar = "Model built with a lasso approach within SAMBA, and a clustering step.\nThe Model is initialized with the true one.",
                                    noCov0="Model built with a lasso approach within SAMBA, and a clustering step.\nThe Model is build whithout excluding covariates.",
                                    standard ="Model built with a lasso approach within SAMBA, and a clustering step."))


  initFolder=paste0("~/Travail/Présentation/Plot/Results")
  if(buildOption!="standard"){
    initFolder=paste0(Folder,"_",opt)
  }
  if(!dir.exists(initFolder)){dir.create(initFolder)}


  initFolder=paste0(initFolder,"/Rsmlx",Rsmlx)
  if(!dir.exists(initFolder)){dir.create(initFolder)}



  # Graphs summarizing building results
  for(opt in buildOption){
    for(meth in buildMethod){
      Folder=paste0(initFolder,"/",meth)
      if(!dir.exists(Folder)){dir.create(Folder)}
      subtitle = paste0(generalsubtitle,Titlelist[[meth]][[opt]])
      for(sim in covariateSize){
        graphsTotalNB(Folder,subtitle,sim,meth,opt,Rsmlx,JPEG,PNG)
        graphsParNB(Folder,subtitle,sim,meth,opt,Rsmlx,JPEG,PNG)
        tableStats(Folder,subtitle,sim,meth,opt,Rsmlx,JPEG,PNG)
      }
    }
  }

  # Comparison graphs
  if(length(buildMethod)>1){
    for(opt in buildOption){
      for(sim in covariateSize){
        Folder=initFolder
        subtitle = generalsubtitle
        graphsCompMethod(Folder,subtitle,sim,buildMethod,opt,Rsmlx,JPEG,compTime,PNG)
        graphsCompMethod2(Folder,subtitle,sim,buildMethod,opt,Rsmlx,JPEG,PNG)
        tableStatsComp(Folder,subtitle,sim,buildMethod,opt,Rsmlx,JPEG,PNG)
      }
    }
  }
}

###############################################################################
tableStats <- function(Folder,subtitle,covariateSize,buildMethod,buildOption="standard",Rsmlx,JPEG,PNG){
  load(paste0("Save/BuildResults_",Rsmlx,".RData"))

  resultModelCov <- resultModel[resultModel$ProjectNumber == paste(covariateSize,"covariates") & resultModel$Method==buildMethod, ]

  errorStatsCov <- errorStats[errorStats$ProjectNumber == paste(covariateSize,"covariates") & errorStats$Method==buildMethod,]

  # Data

  df = data.frame(TypeOfSim=c("cov","corcov"),
                  FDR_mean = c(mean(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"]),
                               mean(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"])),
                  FDR_median = c(median(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"]),
                               median(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"])),
                  FDR_sd = c(sd(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"]),
                             sd(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"])),
                  FDR_q975 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],0.975),
                              quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],0.975)),
                  FDR_q25 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],0.025),
                              quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],0.025)),
                  FNR_mean = c(mean(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"]),
                               mean(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"])),
                  FNR_median = c(median(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"]),
                               median(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"])),
                  FNR_sd = c(sd(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"]),
                             sd(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"])),
                  FNR_q975 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],0.975),
                              quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],0.975)),
                  FNR_q25 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],0.025),
                             quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],0.025)),
                  FN_mean = c(mean(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FN"]),
                               mean(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FN"])),
                  FN_median = c(median(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FN"]),
                                 median(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FN"])),
                  FN_sd = c(sd(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FN"]),
                             sd(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FN"])),
                  FN_q975 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FN"],0.975),
                              quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FN"],0.975)),
                  FN_q25 = c(quantile(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FN"],0.025),
                             quantile(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FN"],0.025)))

  percent <- function(x, digits = 1, format = "f", ...) {
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
  }

  CB <- function(df,char){
    lb = df[,paste0(char,"_q25")]
    ub = df[,paste0(char,"_q975")]
    return(paste0("[",sapply(lb,FUN=function(x){percent(max(0,x))}),
                  ";",sapply(ub,FUN=function(x){percent(min(1,x))}),"]"))}

  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"),FN_CB=CB(df,"FN"))

  tableCov = data.frame(c("With uncorrelated covariates :", "False Discovery Rate","False Negative Rate",
                          paste0("\t • Final Final model without any False Negatives : ",percent(resultModelCov[resultModelCov$TypeOfSim=="cov","NoFNModel"],digits=0)),
                          paste0("\t • Final model is the true one : ",percent(resultModelCov[resultModelCov$TypeOfSim=="cov","TrueModel"],digits=0))),
                        median=c("",percent((as.numeric(format(c(df[df$TypeOfSim=="cov","FDR_median"],df[df$TypeOfSim=="cov","FNR_median"])), scientific=TRUE, digits=2))),"",""),
                        CB=c("",df[df$TypeOfSim=="cov","FDR_CB"],df[df$TypeOfSim=="cov","FNR_CB"],"",""))
  colnames(tableCov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  tableCorcov =data.frame(c("With correlated covariates :","False Discovery Rate","False Negative Rate",
                            paste0("\t • Final Final model without any False Negatives : ",percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"],digits=0)),
                            paste0("\t • Final model is the true one : ",percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"],digits=0))),
                          median=c("",percent(as.numeric(format(c(df[df$TypeOfSim=="corcov","FDR_median"],df[df$TypeOfSim=="corcov","FNR_median"]), scientific=TRUE, digits=2))),"",""),
                          CB=c("",df[df$TypeOfSim=="corcov","FDR_CB"],df[df$TypeOfSim=="corcov","FNR_CB"],"",""))
  colnames(tableCorcov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  table <-  tibble::as_tibble(rbind(tableCov,tableCorcov))




  ft <- flextable(table) %>%
    merge_at(i=1,j=1:3) %>%
    merge_at(i=4,j=1:3) %>%
    merge_at(i=5,j=1:3) %>%
    merge_at(i=6,j=1:3) %>%
    merge_at(i=9,j=1:3) %>%
    merge_at(i=10,j=1:3) %>%
    color(i=c(1,6),color="indianred") %>%
    bold(i=c(1,6)) %>%
    set_table_properties(layout="autofit",width=1) %>%
    bg(i=c(4:5,9:10),j=1:3,bg="#e4e6eb",part="body") %>%
    hline(i = 5, part = "body", border = fp_border_default(color = "grey1", width = 1) ) %>%
    hline(i = 1, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    hline(i = 6, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    vline(j=c(2,1),i=c(2,3,7,8),part="body", border = fp_border_default(color = "grey", width = 1)) %>%
    add_header_lines("Error Rate Comparison Table") %>%
    bold(part="header") %>%
    fontsize(size=22,part="header",i=1)%>%
    fontsize(size=20,part="header",i=2)%>%
    fontsize(size=18,part="body") %>%
    fontsize(size=20,i=c(1,6),part="body") %>%
    align(align = "center", part = "header",i=1) %>%
    fontsize(size=16,i=c(4,5,9,10),part="body") %>%
    add_footer_lines(subtitle) %>%
    fontsize(size=16,part="footer") %>%
    align(align="right",part="footer")

  save_as_html(ft, path = paste0(Folder,"/ErrorTable",covariateSize,".html"),expand=10)
  if(PNG){
    webshot(paste0(Folder,"/ErrorTable",covariateSize,".html"), paste0(Folder,"/ErrorTable",covariateSize,".png"))
  }
  if(JPEG){
    webshot(paste0(Folder,"/ErrorTable",covariateSize,".html"), paste0(Folder,"/ErrorTable",covariateSize,".jpeg"))
  }
  unlink(paste0(Folder,"/ErrorTable",covariateSize,".html"))
}

#########################################################################################

tableStatsComp <- function(Folder,subtitle,covariateSize,buildMethod,buildOption="standard",Rsmlx,JPEG,PNG){
  load(paste0("Save/BuildResults_",Rsmlx,".RData"))

  resultModelCov <- resultModel[resultModel$ProjectNumber == paste(covariateSize,"covariates") & resultModel$Method %in% buildMethod, ]

  errorStatsCov <- errorStats[errorStats$ProjectNumber == paste(covariateSize,"covariates") & errorStats$Method %in% buildMethod,]

# Data
  df = data.frame(TypeOfSim=rep(c("cov","corcov"),each=length(buildMethod)),
                  Method = rep(buildMethod,2),
                  FDR_mean = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = mean)[buildMethod],
                               sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = mean)[buildMethod]),
                  FDR_median = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = median)[buildMethod],
                                 sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = median)[buildMethod]),
                  FDR_sd = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = sd)[buildMethod],
                             sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = sd)[buildMethod]),
                  FDR_q975 = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = function(x){return(quantile(x,0.975))})[paste0(buildMethod,".97.5%")],
                              sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")]),
                  FDR_q25 = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")],
                             sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FDR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")]),
                  FNR_mean = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = mean)[buildMethod],
                               sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = mean)[buildMethod]),
                  FNR_median = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = median)[buildMethod],
                                 sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = median)[buildMethod]),
                  FNR_sd = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = sd)[buildMethod],
                             sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = sd)[buildMethod]),
                  FNR_q975 =  c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")],
                               sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")]),
                  FNR_q25 = c(sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="cov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="cov","Method"]),FUN = function(x){quantile(x,0.05)})[paste0(buildMethod,".2.5%")],
                           sapply(split(errorStatsCov[errorStatsCov$TypeOfSim=="corcov","FNR"],errorStatsCov[errorStatsCov$TypeOfSim=="corcov","Method"]),FUN = function(x){quantile(x,0.05)})[paste0(buildMethod,".2.5%")]))

  percent <- function(x, digits = 1, format = "f", ...) {      # Create user-defined function
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
  }

  CB <- function(df,char){
    lb = df[,paste0(char,"_q25")]
    ub = df[,paste0(char,"_q975")]
    return(paste0("[",sapply(lb,FUN=function(x){percent(max(0,x))}),
                  ";",sapply(ub,FUN=function(x){percent(min(1,x))}),"]"))}

  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"))
  methname=c(lassoSS="Lasso",regression="stepAIC")

  tableCov = data.frame(Rate = c("With uncorrelated covariates :",
                                 paste0("False Discovery Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0("False Negative Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0(" • Final Final model without any False Negatives :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelCov[resultModelCov$TypeOfSim=="cov","NoFNModel"],resultModelCov[resultModelCov$TypeOfSim=="cov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n")),
                                 paste0("  • Final model is the true one :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelCov[resultModelCov$TypeOfSim=="cov","TrueModel"],resultModelCov[resultModelCov$TypeOfSim=="cov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n"))),
                        Median = c("",paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="cov",c("FDR_median")],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],percent),collapse="\n")),
                                   paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="cov",c("FNR_median")],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],percent),collapse="\n")),"",""),
                        CB = c("",paste0("\n",paste0(split(df[df$TypeOfSim=="cov","FDR_CB"],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],collapse="\n")),
                               paste0("\n",paste0(split(df[df$TypeOfSim=="cov","FNR_CB"],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],collapse="\n")),"",""))
  colnames(tableCov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  tableCorcov = data.frame(Rate = c("With correlated covariates :",
                                    paste0("False Discovery Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                    paste0("False Negative Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                    paste0("• Final Final model without any False Negatives :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"],resultModelCov[resultModelCov$TypeOfSim=="corcov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n")),
                                    paste0("• Final model is the true one :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"],resultModelCov[resultModelCov$TypeOfSim=="corcov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n"))),
                           Median = c("",paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="corcov",c("FDR_median")],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],percent),collapse="\n")),
                                      paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="corcov",c("FNR_median")],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],percent),collapse="\n")),"",""),
                           CB = c("",paste0("\n",paste0(split(df[df$TypeOfSim=="corcov","FDR_CB"],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],collapse="\n")),
                                  paste0("\n",paste0(split(df[df$TypeOfSim=="corcov","FNR_CB"],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],collapse="\n")),"",""))
  colnames(tableCorcov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")


  table <-  tibble::as_tibble(rbind(tableCov,tableCorcov))




  ft <- flextable(table) %>%
    merge_at(i=1,j=1:3) %>%
    merge_at(i=4,j=1:3) %>%
    merge_at(i=5,j=1:3) %>%
    merge_at(i=6,j=1:3) %>%
    merge_at(i=9,j=1:3) %>%
    merge_at(i=10,j=1:3) %>%
    color(i=c(1,6),color="indianred") %>%
    bold(i=c(1,6)) %>%
    set_table_properties(layout="autofit",width=1) %>%
    bg(i=c(4:5,9:10),j=1:3,bg="#e4e6eb",part="body") %>%
    hline(i = 5, part = "body", border = fp_border_default(color = "grey1", width = 1) ) %>%
    hline(i = 1, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    hline(i = 6, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    vline(j=c(2,1),i=c(2,3,7,8),part="body", border = fp_border_default(color = "grey", width = 1)) %>%
    add_header_lines("Error Rate Comparison Table") %>%
    bold(part="header") %>%
    fontsize(size=22,part="header",i=1)%>%
    fontsize(size=20,part="header",i=2)%>%
    fontsize(size=18,part="body") %>%
    fontsize(size=20,i=c(1,6),part="body") %>%
    align(align = "center", part = "header",i=1) %>%
    fontsize(size=16,i=c(4,5,9,10),part="body") %>%
    add_footer_lines(subtitle) %>%
    fontsize(size=16,part="footer") %>%
    align(align="right",part="footer")

  save_as_html(ft, path = paste0(Folder,"/ErrorTable",covariateSize,".html"),expand=10)
  webshot(paste0(Folder,"/ErrorTable",covariateSize,".html"), paste0(Folder,"/ErrorTable",covariateSize,".png"))
  unlink(paste0(Folder,"/ErrorTable",covariateSize,".html"))
}


graphsTotalNB <- function(Folder,subtitle,covariateSize,buildMethod,buildOption="standard",Rsmlx,JPEG,PNG){
  load(paste0("Save/BuildResults_",Rsmlx,".RData"))

  fill.vec = c(c("#468b97","#ef6262", "#74C385"),rep("#888888",300))
  gr = "#888888"

  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = paste0(covariateSize," covariates")


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

  covType = paste0("cov")
  corcovType = paste0("corcov")


  # Total number of selection in distinct model :
  valueDisplay = list(cov=c(cAGE = resultCovariateCov[resultCovariateCov$Covariate=="cAGE" &
                                                        resultCovariateCov$TypeOfSim==covType,"SelectedInDistinctModel"]*100,
                            RACE = resultCovariateCov[resultCovariateCov$Covariate=="RACE" &
                                                        resultCovariateCov$TypeOfSim==covType,"SelectedInDistinctModel"]*100,
                            SEX = resultCovariateCov[resultCovariateCov$Covariate=="SEX" &
                                                       resultCovariateCov$TypeOfSim==covType,"SelectedInDistinctModel"]*100,
                            FP = round(mean(resultCovariateCov[!(resultCovariateCov$Covariate %in% c("cAGE","RACE","SEX"))
                                                               & resultCovariateCov$TypeOfSim==covType,"SelectedInDistinctModel"]),digits = 3)*100),
                      corcov=c(cAGE = resultCovariateCov[resultCovariateCov$Covariate=="cAGE" &
                                                           resultCovariateCov$TypeOfSim==corcovType,"SelectedInDistinctModel"]*100,
                               RACE = resultCovariateCov[resultCovariateCov$Covariate=="RACE" &
                                                           resultCovariateCov$TypeOfSim==corcovType,"SelectedInDistinctModel"]*100,
                               SEX = resultCovariateCov[resultCovariateCov$Covariate=="SEX" &
                                                          resultCovariateCov$TypeOfSim==corcovType,"SelectedInDistinctModel"]*100,
                               FP = round(mean(resultCovariateCov[!(resultCovariateCov$Covariate %in% c("cAGE","RACE","SEX"))
                                                                  & resultCovariateCov$TypeOfSim==corcovType,"SelectedInDistinctModel"]),digits = 3)*100))



  cov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,],aes(x=factor(Covariate,levels=orderList$cov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=c(rep(0.5,3),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"])) - 3)),
             color=fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"]))],
             fill = fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"]))])+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With uncorrelated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim==covType,"NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim==covType,"TrueModel"]*100,"%"))+
    geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
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
    annotate("text",x=4,y=0.9,hjust=0,color="#468b97",size=5,
             label=paste0("cAGE :",valueDisplay$cov[["cAGE"]],"%"))+
    annotate("text",x=4,y=0.75,hjust=0,color="#ef6262",size=5,
             label=paste0("RACE :",valueDisplay$cov[["RACE"]],"%"))+
    annotate("text",x=4,y=0.6,hjust=0,color="#74C385",size=5,
             label=paste0("SEX :",valueDisplay$cov[["SEX"]],"%"))+
    annotate("text",x=covariateSize,y=valueDisplay$cov[["FP"]]/100+0.1,hjust=1,color="#9E9FA5",fontface = 'italic',size=5,label=paste0(valueDisplay$cov["FP"],"%"))+
    geom_segment(x=4,y=valueDisplay$cov[["FP"]]/100,xend=covariateSize+0.5,yend=valueDisplay$cov[["FP"]]/100,color="#9E9FA5",linewidth=1,linetype="dashed")+
    coord_cartesian(ylim=c(0,1),clip="off")

  corcov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,],aes(x=factor(Covariate,levels=orderList$corcov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=c(rep(0.5,3),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"]))-3)),
             color=fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"]))],
             fill = fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"]))])+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With correlated covariate, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim==corcovType,"NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim==corcovType,"TrueModel"]*100,"%"))+
    geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
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
    annotate("text",x=4,y=0.9,hjust=0,color="#468b97",size=5,
             label=paste0("cAGE :",valueDisplay$corcov[["cAGE"]],"%"))+
    annotate("text",x=4,y=0.75,hjust=0,color="#ef6262",size=5,
             label=paste0("RACE :",valueDisplay$corcov[["RACE"]],"%"))+
    annotate("text",x=4,y=0.6,hjust=0,color="#74C385",size=5,
             label=paste0("SEX :",valueDisplay$corcov[["SEX"]],"%"))+
    annotate("text",x=covariateSize,y=valueDisplay$corcov[["FP"]]/100+0.15,hjust=1,color="#9E9FA5",fontface = 'italic',size=5,label=paste0(valueDisplay$corcov["FP"],"%"))+
    geom_segment(x=4,y=valueDisplay$corcov[["FP"]]/100,xend=covariateSize+.5,yend=valueDisplay$corcov[["FP"]]/100,color="#9E9FA5",linewidth=1,linetype="dashed")+
    coord_cartesian(ylim=c(0,1),clip="off")

  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )

  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,".png"),
           height = 2500, width = 5500+covariateSize, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,".jpeg"),
           height = 2500, width = 5500+covariateSize, units = "px",device=grDevices::jpeg)
  }


  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"Cov.png"),
           height = 1500, width = 5500+covariateSize, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"Cov.jpeg"),
           height = 1500, width = 5500+covariateSize, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"Corcov.png"),
           height = 1500, width = 5500+covariateSize, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"Corcov.jpeg"),
           height = 1500, width = 5500+covariateSize, units = "px",device=grDevices::jpeg)
  }

}

###########################################################################################################################################

graphsParNB <- function(Folder,subtitle,covariateSize, buildMethod,buildOption,Rsmlx,JPEG,PNG){
  load(paste0("Save/BuildParResults_",Rsmlx,".RData"))

  fill.vec = c(c("#468b97","#ef6262", "#74C385"),rep("#888888",200))
  gr = "#888888"

  covType = paste0("cov")
  corcovType = paste0("corcov")

  name = list("varphi[S]" = latex2exp::TeX(r"($\varphi_S$)",output="character"),
           "varphi[L]" = latex2exp::TeX(r"($\varphi_L$)",output="character"),
           "delta[Ab]" = latex2exp::TeX(r"($\delta_{Ab}$)",output="character"))


  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = latex2exp::TeX(paste0(covariateSize," covariates"),output="character")


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


  covariateAB =list(cov =  unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==covType
                                                             & CovariateModelSelectionCov$Parameter=="delta[Ab]","Covariate"]),
                    corcov = unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==corcovType
                                                               & CovariateModelSelectionCov$Parameter=="delta[Ab]","Covariate"]))

  covariateL = list(cov =  unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==covType
                                                             & CovariateModelSelectionCov$Parameter=="varphi[L]","Covariate"]),
                    corcov = unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==corcovType
                                                               & CovariateModelSelectionCov$Parameter=="varphi[L]","Covariate"]))

  covariateS = list(cov =  unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==covType
                                                             & CovariateModelSelectionCov$Parameter=="varphi[S]","Covariate"]),
                    corcov = unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==corcovType
                                                               & CovariateModelSelectionCov$Parameter=="varphi[S]","Covariate"]))
  lim = list(cov=length(union(covariateAB$cov,union(covariateL$cov,covariateS$cov))),
             corcov = length(union(covariateAB$corcov,union(covariateL$corcov,covariateS$corcov))))


  repAB = list(cov=2,corcov=2)
  repL=list(cov=list(av=1,ap=1),corcov=list(av=1,ap=1))
  repS = list(cov=2,corcov=2)
  if(!("cAGE" %in% covariateAB$cov) && !("RACE" %in% covariateAB$cov)){
    repAB$cov = 0
  }else if(!("cAGE" %in% covariateAB$cov) || !("RACE" %in% covariateAB$cov)){
    repAB$cov = 1
  }
  if(!("cAGE" %in% covariateL$cov)){
    repL$cov$av = 0
  }
  if(!("SEX" %in% covariateL$cov)){
    repL$cov$ap = 0
  }
  if(!("SEX" %in% covariateS$cov) && !("RACE" %in% covariateS$cov)){
    repS$cov = 0
  }else if(!("SEX" %in% covariateS$cov) || !("RACE" %in% covariateS$cov)){
    repS$cov = 1
  }



  if(!("cAGE" %in% covariateAB$corcov) && !("RACE" %in% covariateAB$corcov)){
    repAB$corcov = 0
  }else if(!("cAGE" %in% covariateAB$corcov) || !("RACE" %in% covariateAB$corcov)){
    repAB$corcov = 1
  }
  if(!("cAGE" %in% covariateL$corcov)){
    repL$corcov$av = 0
  }
  if(!("SEX" %in% covariateL$corcov)){
    repL$corcov$ap = 0
  }
  if(!("SEX" %in% covariateS$corcov) && !("RACE" %in% covariateS$corcov)){
    repS$corcov = 0
  }else if(!("SEX" %in% covariateS$corcov) || !("RACE" %in% covariateS$corcov)){
    repS$corcov = 1
  }

  ## Selection In Model

  valueDisplay = data.frame()
  for(t in c(covType,corcovType)){
      value=c(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="cAGE","ProportionSelected"]*100,
              resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="RACE","ProportionSelected"]*100,
              resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="SEX","ProportionSelected"]*100,
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="cAGE","ProportionSelected"]*100),digits=2),
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="RACE","ProportionSelected"]*100),digits=2),
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="SEX","ProportionSelected"]*100),digits=2))
      valuemax = c(max(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="cAGE","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm = TRUE),
                   max(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="RACE","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm = TRUE),
                   max(resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="SEX","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm = TRUE))

      valueDisplay <- rbind(valueDisplay,data.frame(cov=c(rep("cov",3),rep("FP",3)),
                                                    value=value,
                                                    type=t,
                                                    coor=c(2,3,4,valuemax),
                                                    Parameter=rep(c("varphi[S]","varphi[L]","delta[Ab]"),2),
                                                    distext=paste0(c("cAGE :","RACE :","SEX :",rep("",3)),value,"%")))
  }


  cov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==covType,],aes(x=factor(Covariate,levels=orderList$cov)))+
    geom_bar(position=position_dodge(preserve = "single"),
             alpha=c(c(rep(0.2,repAB$cov),0.5,
                       rep(0.2,length(covariateAB$cov) - 1 - repAB$cov)),

                     c(rep(0.2,repL$cov$av),0.5,rep(0.2,repL$cov$ap),
                       rep(0.2,length(covariateL$cov)- 1 - repL$cov$av - repL$cov$ap )),

                     c(0.5,rep(0.2,repS$cov),
                       rep(0.2,length(covariateS$cov)- 1 -repS$cov))),

             color=c(c(rep(gr,repAB$cov),fill.vec[3],
                       rep(gr,length(covariateAB$cov)- 1 - repAB$cov)),

                     c(rep(gr,repL$cov$av),fill.vec[2],rep(gr,repL$cov$ap),
                       rep(gr,length(covariateL$cov)-1 - repL$cov$av - repL$cov$ap)),

                     c(fill.vec[1],rep(gr,repS$cov),
                       rep(gr,length(covariateS$cov)-1 -repS$cov))),


             fill = c(c(rep(gr,repAB$cov),fill.vec[3],
                        rep(gr,length(covariateAB$cov)- 1 - repAB$cov)),

                      c(rep(gr,repL$cov$av),fill.vec[2],rep(gr,repL$cov$ap),
                        rep(gr,length(covariateL$cov)-1 - repL$cov$av - repL$cov$ap)),

                      c(fill.vec[1],rep(gr,repS$cov),
                        rep(gr,length(covariateS$cov)-1 -repS$cov))))+
    facet_grid(Parameter~.,labeller=label_parsed)+
    xlab("Covariate")+
    ylab("Count")+
    ggtitle("With uncorrelated covariates, ",
            subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim==covType,"NoFNModel"]*100,"%","\n",
                             "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim==covType,"TrueModel"]*100,"%"))+
    geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    ylim(c(0,100))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=90,hjust=0),color=c("#468b97","#ef6262", "#74C385"),size=5)+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$cov,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
  coord_cartesian(ylim=c(0,100),clip="off")

  corcov =  ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==corcovType,],aes(x=factor(Covariate,levels=orderList$corcov)))+
      geom_bar(position=position_dodge(preserve = "single"),
               alpha=c(c(rep(0.2,repAB$corcov),0.5,
                         rep(0.2,length(covariateAB$corcov) - 1 - repAB$corcov)),

                       c(rep(0.2,repL$corcov$av),0.5,rep(0.2,repL$corcov$ap),
                         rep(0.2,length(covariateL$corcov)- 1 - repL$corcov$av - repL$corcov$ap )),

                       c(0.5,rep(0.2,repS$corcov),
                         rep(0.2,length(covariateS$corcov)- 1 -repS$corcov))),

               color=c(c(rep(gr,repAB$corcov),fill.vec[3],
                         rep(gr,length(covariateAB$corcov)- 1 - repAB$corcov)),

                       c(rep(gr,repL$corcov$av),fill.vec[2],rep(gr,repL$corcov$ap),
                         rep(gr,length(covariateL$corcov)-1 - repL$corcov$av - repL$corcov$ap)),

                       c(fill.vec[1],rep(gr,repS$corcov),
                         rep(gr,length(covariateS$corcov)-1 -repS$corcov))),


               fill = c(c(rep(gr,repAB$corcov),fill.vec[3],
                          rep(gr,length(covariateAB$corcov)- 1 - repAB$corcov)),

                        c(rep(gr,repL$corcov$av),fill.vec[2],rep(gr,repL$corcov$ap),
                          rep(gr,length(covariateL$corcov)-1 - repL$corcov$av - repL$corcov$ap)),

                        c(fill.vec[1],rep(gr,repS$corcov),
                          rep(gr,length(covariateS$corcov)-1 -repS$corcov))))+
      facet_grid(Parameter~.,labeller=label_parsed)+
      xlab("Covariate")+
      ylab("Count")+
      ggtitle("With correlated covariates, ",
              subtitle= paste0("   • Final Final model without any False Negatives : ",resultModelCov[resultModelCov$TypeOfSim==corcovType,"NoFNModel"]*100,"%","\n",
                               "   • Final model is the true one  : ",resultModelCov[resultModelCov$TypeOfSim==corcovType,"TrueModel"]*100,"%"))+
      geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 6, angle = 90))+
      ylim(c(0,100))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=14))+
      theme(strip.text = element_text(size = 16))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.text = element_text(size=14))+
      theme(legend.title = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))+
      theme(plot.subtitle = element_text(size=12))+
      geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=90,hjust=0),color=c("#468b97","#ef6262", "#74C385"),size=5)+
      geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$corcov,y=coor+10,vjust=0,hjust=1),color="#9E9FA5",fontface = 'italic',size=5)+
      coord_cartesian(ylim=c(0,100),clip="off")

  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )

  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,".png"),
           height = 2500, width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,".jpeg"),
           height = 2500, width = 5500, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"Cov.png"),
           height = 1500, width = 5500+covariateSize, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"Cov.jpeg"),
           height = 1500, width = 5500+covariateSize, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )
  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"Corcov.png"),
           height = 1500, width = 5500+covariateSize, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"Corcov.jpeg"),
           height = 1500, width = 5500+covariateSize, units = "px",device=grDevices::jpeg)
  }


}
#############################################################################################################
graphsCompMethod <- function(Folder,subtitle,covariateSize,buildMethod,buildOption,Rsmlx,JPEG,compTime=TRUE,PNG){
  library(ggplot2,quietly=TRUE)
  source("~/Travail/00_Theme.R")

  fill.vec = c(c("#468b97","#ef6262", "#74C385"),rep("#888888",200))
  gr = "#888888"


  load(paste0("Save/BuildResults_",Rsmlx,".RData"))


  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = paste0(covariateSize," covariates")

  covType = paste0("cov")
  corcovType = paste0("corcov")

  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method %in% buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber == covariateSizeCar,]
  resultCovariateCov <- resultCovariate[resultCovariate$Method %in% buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]
  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method %in% buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]

  computationStatsCov <- computationStats[computationStats$Method %in% buildMethod
                                          & computationStats$ProjectNumber==covariateSizeCar,]

  resultModelCov <- resultModel[resultModel$Method %in% buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]
  resultModelCov <- cbind(resultModelCov,text=paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))

  errorStatsCov <- errorStats[errorStats$Method %in% buildMethod,]

  # Total number of selection in distinct model :
  lim=data.frame()
  for(t in c(covType,corcovType)){
    for(meth in buildMethod){
      lim = rbind(lim,data.frame(TypeOfSim=t,
                                 lim=length(unique(CovariateModelSelectionCov[CovariateModelSelectionCov$Method==meth
                                                                              & CovariateModelSelectionCov$TypeOfSim==t,
                                                                              "Covariate"])),
                                 Method=meth))
    }
  }

  valueDisplay = data.frame()
  for(t in c(covType,corcovType)){
    for(meth in buildMethod){
      value=c(resultCovariateCov[resultCovariateCov$Covariate=="cAGE"
                                 & resultCovariateCov$Method == meth
                                 & resultCovariateCov$TypeOfSim==t,"SelectedInDistinctModel"]*100,
              resultCovariateCov[resultCovariateCov$Covariate=="RACE"
                                 & resultCovariateCov$Method == meth
                                 & resultCovariateCov$TypeOfSim==t,"SelectedInDistinctModel"]*100,

              resultCovariateCov[resultCovariateCov$Covariate=="SEX"
                                 & resultCovariateCov$Method == meth
                                 & resultCovariateCov$TypeOfSim==t,"SelectedInDistinctModel"]*100,
              round(mean(resultCovariateCov[!(resultCovariateCov$Covariate %in% c("cAGE","RACE","SEX"))
                                            & resultCovariateCov$Method == meth
                                            & resultCovariateCov$TypeOfSim==t,"SelectedInDistinctModel"]),digits = 3)*100)

      valuemax = max(resultCovariateCov[!(resultCovariateCov$Covariate %in% c("cAGE","RACE","SEX"))
                                            & resultCovariateCov$Method == meth
                                            & resultCovariateCov$TypeOfSim==t,"SelectedInDistinctModel"][(lim[lim$TypeOfSim==t & lim$Method==meth,"lim"]-5):lim[lim$TypeOfSim==t & lim$Method==meth,"lim"]]*100,na.rm = TRUE)

      valueDisplay <- rbind(valueDisplay,data.frame(cov=c(rep("cov",3),"FP"),
                                                    value=value,
                                                    Method=meth,
                                                    type=t,
                                                    coor=valuemax,
                                                    Parameter=c("varphi[S]","varphi[L]","delta[Ab]","NONE"),
                                                    distext=paste0(c("cAGE :","RACE :","SEX :",""),value,"%")))
    }
  }
  valueModel = cbind(resultModelCov,Parameter=c("delta[Ab]"),text = paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))


  cov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,],aes(x=factor(Covariate,levels=orderList$cov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=rep(c(rep(0.5,3),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"])) - 3)),length(buildMethod)),
             color=rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"]))],length(buildMethod)),
             fill =rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==covType,"Covariate"]))],length(buildMethod)))+
    facet_grid(factor(Method,levels=c("regression","lassoSS","lassoSSCl","ClustOfVar"))~.,labeller = as_labeller(c(lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step.",regression="StepAIC")))+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With uncorrelated covariates, ")+
    geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    theme(strip.text.x = element_text(size = 8))+
    ylim(c(0,1))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 14))+
    theme(plot.title = element_text(size=24,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=16))+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="varphi[S]",],
              mapping=aes(label=distext,x=4,y=0.9),color="#468b97",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="varphi[L]",],
              mapping=aes(label=distext,x=4,y=0.75),color="#ef6262",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="delta[Ab]",],
              mapping=aes(label=distext,x=4,y=0.6),color="#74C385",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=max(lim[lim$TypeOfSim==covType,"lim"]),y=(coor+10)/100,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    geom_text(data = resultModelCov[resultModelCov$TypeOfSim==covType,],mapping=aes(label=text),x=covariateSize,y=1,hjust=1,size=5,vjust=1)+
    coord_cartesian(ylim=c(0,1),clip="off")

  corcov = ggplot(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,],aes(x=factor(Covariate,levels=orderList$corcov), y = SelectedInDistinctModel))+
    geom_bar(stat="identity",position="dodge",
             alpha=rep(c(rep(0.5,3),rep(0.2,length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"])) - 3)),length(buildMethod)),
             color=rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"]))],length(buildMethod)),
             fill =rep(fill.vec[1:length(unique(resultCovariateCov[resultCovariateCov$TypeOfSim==corcovType,"Covariate"]))],length(buildMethod)))+
    facet_grid(factor(Method,levels=c("regression","lassoSS","lassoSSCl","ClustOfVar"))~.,labeller = as_labeller(c(lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step.",regression="StepAIC")))+
    xlab("Covariate")+
    ylab("Frequency")+
    ggtitle("With correlated covariates, ")+
    geom_segment(x=3.5,y=0,xend=3.5,yend=1,color="#862B0D",linewidth=1,linetype="twodash")+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 10))+
    theme(strip.text.x = element_text(size = 8))+
    ylim(c(0,1))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 14))+
    theme(plot.title = element_text(size=24,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=16))+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="varphi[S]",],
              mapping=aes(label=distext,x=4,y=0.9),color="#468b97",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="varphi[L]",],
              mapping=aes(label=distext,x=4,y=0.75),color="#ef6262",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="cov" & valueDisplay$Parameter=="delta[Ab]",],
              mapping=aes(label=distext,x=4,y=0.6),color="#74C385",size=5,hjust=0)+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=max(lim[lim$TypeOfSim==corcovType,"lim"]),y=(coor+10)/100,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5)+
    geom_text(data = resultModelCov[resultModelCov$TypeOfSim==corcovType,],mapping=aes(label=text),x=covariateSize,y=1,hjust=1,size=5,vjust=1)+
    coord_cartesian(ylim=c(0,1),clip="off")


  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
                 height = 4500, width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }
  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 4500, width = 5500, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )

  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = 2500, width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = 2500, width = 5500, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = 2500, width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){

    ggsave(paste0(Folder,"/NumberOfSelection",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = 2500, width = 5500, units = "px",device=grDevices::jpeg)
  }

  colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]

  if(compTime){
    ## Computation Stats
    valueDisplayTime = data.frame()
    for(t in c(covType,corcovType)){
      for(meth in buildMethod){
        aux = computationStatsCov[computationStatsCov$TypeOfSim==t & computationStatsCov$Method==meth,]
        valueDisplayTime = rbind(valueDisplayTime,data.frame(Method=meth,
                                                             TypeOfSim = t,
                                                             Mean = median(aux$time),
                                                             textMean = paste0(round(median(aux$time)), " s"),
                                                             textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                             StandardDeviation = sd(aux$time)))
      }
    }
    valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))

    ggplot(computationStatsCov[computationStatsCov$TypeOfSim %in% c(covType,corcovType),],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=time))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=rep(col,2),
                   fill=rep(colpas,2))+
      geom_text(data=valueDisplayTime, mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6,color=rep(colFonce,2))+
      facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates",cov2="With uncorrelated covariates",corcov2="With correlated covariates")))+
      #stat_boxplot(geom = "errorbar") +
      xlab("Method used")+
      ylab("Computation time (s)")+
      ggtitle("Computation Time Comparison",
              subtitle=subtitle)+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))

    if(PNG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
             height = 2800, width = 2500, units = "px", bg='transparent',device=grDevices::png)
    }

    if(JPEG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
             height = 2800, width = 2500, units = "px",device=grDevices::jpeg)
    }


    ggplot(computationStatsCov[computationStatsCov$TypeOfSim ==covType,],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=time))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=col,
                   fill=colpas)+
      geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim==covType,], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6,color=colFonce)+
      #stat_boxplot(geom = "errorbar") +
      xlab("Method used")+
      ylab("Computation time (s)")+
      ggtitle("Computation Time Comparison",
              subtitle=paste0(subtitle,"With uncorrelated covariates,"))+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))

    if(PNG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
             height = 2000, width = 2000, units = "px", bg='transparent',device=grDevices::png)
    }

    if(JPEG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
             height = 2000, width = 2000, units = "px",device=grDevices::jpeg)
    }

    ggplot(computationStatsCov[computationStatsCov$TypeOfSim ==corcovType,],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=time))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=col,
                   fill=colpas)+
      geom_text(data=valueDisplayTime[valueDisplayTime$TypeOfSim==corcovType,], mapping=aes(label=textMean,y=Mean,x=coor),vjust=-0.2,fontface="bold",size=6,color=colFonce)+
      #stat_boxplot(geom = "errorbar") +
      xlab("Method used")+
      ylab("Computation time (s)")+
      ggtitle("Computation Time Comparison",
              subtitle=paste0(subtitle,"With correlated covariates,"))+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x)) +
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))

    if(PNG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
                   height = 2000, width = 2000, units = "px", bg='transparent',device=grDevices::png)
      }

    if(JPEG){
      ggsave(paste0(Folder,"/ComputationTime",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
             height = 2000, width = 2000, units = "px",device=grDevices::jpeg)
    }





    ## Iteration time
    valueDisplayTime = data.frame()
    for(t in c(covType,corcovType)){
      for(meth in buildMethod){
        aux = computationStatsCov[computationStatsCov$TypeOfSim==t & computationStatsCov$Method==meth,]
        valueDisplayTime = rbind(valueDisplayTime,data.frame(Method=meth,
                                                             TypeOfSim = t,
                                                             Mean = median(aux$iteration),
                                                             textMean = paste0(median(aux$iteration)),
                                                             textStandardDeviation = paste0("sd = ", round(sd(aux$iteration),digits=2)),
                                                             StandardDeviation = sd(aux$iteration)))
      }
    }

    valueDisplayTime <- cbind(valueDisplayTime,coor=1:length(buildMethod))

    ggplot(computationStatsCov[computationStatsCov$TypeOfSim %in% c(covType,corcovType),],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=iteration))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=rep(col,2),
                   fill=rep(colpas,2))+
      facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates",cov2="With uncorrelated covariates",corcov2="With correlated covariates")))+
      xlab("Method used")+
      ylab("Iteration Count")+
      ggtitle("Number of Iterations Comparison",
              subtitle=subtitle)+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))


    if(PNG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
             height = 2800, width = 2500, units = "px", bg='transparent',device=grDevices::png)
    }


    if(JPEG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
             height = 2800, width = 2500, units = "px",device=grDevices::jpeg)
    }

    ggplot(computationStatsCov[computationStatsCov$TypeOfSim == covType,],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=iteration))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=col,
                   fill=colpas)+
      xlab("Method used")+
      ylab("Iteration Count")+
      ggtitle("Number of Iterations Comparison",
              subtitle=paste0(subtitle,"With uncorrelated covariates,"))+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))


    if(PNG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
             height = 2000, width = 2000, units = "px", bg='transparent',device=grDevices::png)
    }


    if(JPEG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
             height = 2000, width = 2000, units = "px",device=grDevices::jpeg)
    }

    ggplot(computationStatsCov[computationStatsCov$TypeOfSim == corcovType,],aes(x=factor(Method,levels=c("regression","lassoSS","lassoSSCl")),y=iteration))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=col,
                   fill=colpas)+
      xlab("Method used")+
      ylab("Iteration Count")+
      ggtitle("Number of Iterations Comparison",
              subtitle=paste0(subtitle,"With correlated covariates,"))+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSSCl="Lasso with\nclustering step",lassoSS="Lasso",ClustOfVar="ClustOfVar"))+
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))


    if(PNG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
             height = 2000, width = 2000, units = "px", bg='transparent',device=grDevices::png)
    }

    if(JPEG){
      ggsave(paste0(Folder,"/IterationCount",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
             height = 2000, width = 2000, units = "px",device=grDevices::jpeg)
    }
  }


  colFonce = colFonce = c("#5c6e39","#563f61","#703527","#024154","#524b43")[1:length(buildMethod)]
  col = c("#a6c46a","#8e6aa0","#ee6c4d","#007194","#9D8F80","#FFD447")[1:length(buildMethod)]
  colpas = c("#e0e6c6","#d0c1d7","#f8c2b4","#99e7ff","#cac2ba","#ffe591")[1:length(buildMethod)]## FDR

  valueDisplayFDR = data.frame()
  for(t in c(covType,corcovType)){
    for(meth in buildMethod){
      aux = errorStatsCov[errorStatsCov$TypeOfSim==t & errorStatsCov$Method==meth,]
      valueDisplayFDR = rbind(valueDisplayFDR,data.frame(Method=meth,
                                                          TypeOfSim = t,
                                                          Mean = mean(aux$FDR),
                                                          Median = median(aux$FDR),
                                                         text = paste0("Mean : ",round(mean(aux$FDR),digits=2),"\nMedian : ",round(median(aux$FDR),digits=2)),
                                                          textStandardDeviation = paste0("sd = ", round(sd(aux$time),digits=2)),
                                                           StandardDeviation = sd(aux$time)))
    }
  }

  valueDisplayFDR <- cbind(valueDisplayFDR,coor=1:length(buildMethod))

  # FDR =
    ggplot(errorStatsCov[errorStatsCov$TypeOfSim %in% c(covType,corcovType),],aes(y=FDR,x=factor(Method,levels=c("regression","lassoSS","lassoSSCl"))))+
      geom_boxplot(lwd=0.5,alpha=0.6,
                   color=rep(col,2),
                   fill=rep(colpas,2))+
    geom_text(data=valueDisplayFDR, mapping=aes(label=text,x=coor),y=0.10,vjust=-0.2,fontface="bold",size=3,color=rep(colFonce,2))+
    facet_grid(.~TypeOfSim,labeller = as_labeller(c(cov="With uncorrelated covariates",corcov="With correlated covariates",cov2="With uncorrelated covariates",corcov2="With correlated covariates"))) +
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
      xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=subtitle)+
      scale_x_discrete(labels=c(regression="stepAIC",lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step"))+
      scale_fill_manual(values=cbPalette)+
      theme(axis.text.x = element_text(size = 10))+
      theme(axis.text.y = element_text(size = 8))+
      theme(axis.title = element_text(size=12))+
      theme(strip.text = element_text(size = 12))+
      theme(plot.subtitle = element_text(size=12))+
      theme(plot.title = element_text(size=16,color="#ee6c4d"))




  if(PNG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = 1700, width = 1700, units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 1700, width = 1700, units = "px",device=grDevices::jpeg)
  }

  ggplot(errorStatsCov[errorStatsCov$TypeOfSim == covType,],aes(y=FDR,x=factor(Method,levels=c("regression","lassoSS","lassoSSCl"))))+
    geom_boxplot(lwd=0.5,alpha=0.6,
                 color=col,
                 fill=colpas)+
    geom_text(data=valueDisplayFDR[valueDisplayFDR$TypeOfSim==covType,], mapping=aes(label=text,x=coor),y=0.10,vjust=-0.2,fontface="bold",size=4,color=colFonce)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=paste0(subtitle,"With uncorrelated covariates,"))+
    scale_x_discrete(labels=c(regression="stepAIC",lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step"))+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))

  if(PNG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = 1700, width = 1700, units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = 1700, width = 1700, units = "px",device=grDevices::jpeg)
  }

  ggplot(errorStatsCov[errorStatsCov$TypeOfSim == corcovType,],aes(y=FDR,x=factor(Method,levels=c("regression","lassoSS","lassoSSCl"))))+
    geom_boxplot(lwd=0.5,alpha=0.6,
                 color=col,
                 fill=colpas)+
    geom_text(data=valueDisplayFDR[valueDisplayFDR$TypeOfSim==corcovType,], mapping=aes(label=text,x=coor),y=0.10,vjust=-0.2,fontface="bold",size=4,color=colFonce)+
    guides(fill = guide_legend(title = "Method Used :"), color= guide_legend(title = "Method Used :")) +
    xlab("Method used")+
    ggtitle("False Discovery Rate Distribution",
            subtitle=paste0(subtitle,"With correlated covariates,"))+
    scale_x_discrete(labels=c(regression="stepAIC",lassoSS="Lasso",ClustOfVar="ClustOfVar",lassoSSCl="Lasso with\nclustering step"))+
    scale_fill_manual(values=cbPalette)+
    theme(axis.text.x = element_text(size = 10))+
    theme(axis.text.y = element_text(size = 8))+
    theme(axis.title = element_text(size=12))+
    theme(strip.text = element_text(size = 12))+
    theme(plot.subtitle = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))

  if(PNG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = 1700, width = 1700, units = "px", bg='transparent',device=grDevices::png)
  }

  if(JPEG){
    ggsave(paste0(Folder,"/ErrorControl",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = 1700, width = 1700, units = "px",device=grDevices::jpeg)
  }

}

#############################################################################################################
graphsCompMethod2 <- function(Folder,subtitle,covariateSize,buildMethod,buildOption,Rsmlx,JPEG,PNG){
  load(paste0("Save/BuildParResults_",Rsmlx,".RData"))

  fill.vec = c(c("#468b97","#ef6262", "#74C385"),rep("#888888",200))
  gr = "#888888"

  covType = paste0("cov")
  corcovType = paste0("corcov")

  newbuildMethod <- c()
  for(k in 1:length(buildMethod)){
    if(buildMethod[k]=="regression"){
      newbuildMethod[k] <- "StepAIC"
    }else if(buildMethod[k]=="lassoSSCl"){
      newbuildMethod[k] <- "Lasso with\nclustering step."
    }else if(buildMethod[k]=="lassoSS"){
      newbuildMethod[k] <- "Lasso"
    }else if(buildMethod[k]=="CoVSS"){
      newbuildMethod[k] <- "ClustOfVar"
    }
  }


  covToKeep = union(orderList$cov[1:covariateSize],orderList$corcov[1:covariateSize])
  covariateSizeCar = latex2exp::TeX(paste0(covariateSize," covariates"),output="character")


  CovariateModelSelectionCov <- CovariateModelSelection[CovariateModelSelection$Method %in% buildMethod
                                                        & CovariateModelSelection$Covariate %in% covToKeep
                                                        & CovariateModelSelection$ProjectNumber == covariateSizeCar,]

  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="regression","Method"] <- "StepAIC"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="lassoSS","Method"] <- "Lasso"
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="lassoSSCl","Method"] <- "Lasso with\nclustering step."
  CovariateModelSelectionCov[CovariateModelSelectionCov$Method=="CoVSS","Method"] <- "ClustOfVar"


  resultCovariateCov <- resultCovariate[resultCovariate$Method %in% buildMethod
                                        & resultCovariate$Covariate %in% covToKeep
                                        & resultCovariate$ProjectNumber==covariateSizeCar,]

  resultCovariateCov[resultCovariateCov$Method=="regression","Method"] <- "StepAIC"
  resultCovariateCov[resultCovariateCov$Method=="lassoSS","Method"] <- "Lasso"
  resultCovariateCov[resultCovariateCov$Method=="lassoSSCl","Method"] <- "Lasso with\nclustering step."
  resultCovariateCov[resultCovariateCov$Method=="CoVSS","Method"] <- "ClustOfVar"


  resultCovariateParCov <- resultCovariatePar[resultCovariatePar$Method %in% buildMethod
                                              & resultCovariatePar$Covariate %in% covToKeep
                                              & resultCovariatePar$ProjectNumber==covariateSizeCar,]

  resultCovariateParCov[resultCovariateParCov$Method=="regression","Method"] <- "StepAIC"
  resultCovariateParCov[resultCovariateParCov$Method=="lassoSS","Method"] <- "Lasso"
  resultCovariateParCov[resultCovariateParCov$Method=="lassoSSCl","Method"] <- "Lasso with\nclustering step."
  resultCovariateParCov[resultCovariateParCov$Method=="CoVSS","Method"] <- "ClustOfVar"

  resultModelCov <- resultModel[resultModel$Method %in% buildMethod
                                & resultModel$ProjectNumber==covariateSizeCar,]
  resultModelCov <- cbind(resultModelCov,text=paste0("Model without False Negatives : ",resultModelCov$NoFNModel*100,"%"))

  resultModelCov[resultModelCov$Method=="regression","Method"] <- "StepAIC"
  resultModelCov[resultModelCov$Method=="lassoSS","Method"] <- "Lasso"
  resultModelCov[resultModelCov$Method=="lassoSSCl","Method"] <- "Lasso with\nclustering step."
  resultModelCov[resultModelCov$Method=="CoVSS","Method"] <- "ClustOfVar"


  buildMethod <- newbuildMethod

  covariate=list()
  for(t in c(covType,corcovType)){
    aux = list()
    for(meth in buildMethod){
      covmetht = list(S=unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==t
                                                          & CovariateModelSelectionCov$Method==meth
                                                          & CovariateModelSelectionCov$Parameter=="varphi[S]"
                                                          ,"Covariate"]),
                      L=unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==t
                                                          & CovariateModelSelectionCov$Method==meth
                                                          & CovariateModelSelectionCov$Parameter=="varphi[L]"
                                                          ,"Covariate"]),
                      AB=unique(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==t
                                                           & CovariateModelSelectionCov$Method==meth
                                                           & CovariateModelSelectionCov$Parameter=="delta[Ab]"
                                                           ,"Covariate"]))
      aux <- append(aux,list(covmetht))
      names(aux)[length(aux)] <- meth
    }
    covariate <- append(covariate,list(aux))
    names(covariate)[length(covariate)] <- t
  }

  lim = list(cov=length(unique(unlist(covariate$cov))),
             corcov = length(unique(unlist(covariate$cov))))

  rep = data.frame()
  for(t in c(covType,corcovType)){
    for(meth in buildMethod){
      repAB = 2
      if(!("cAGE" %in% covariate[[t]][[meth]]$AB) && !("RACE" %in% covariate[[t]][[meth]]$AB)){
        repAB = 0
      }else if(!("cAGE" %in% covariate[[t]][[meth]]$AB) || !("RACE" %in% covariate[[t]][[meth]]$AB)){
        repAB = 1
      }
      rep <- rbind(rep,data.frame(TypeOfSim=t,Method=meth,par="AB",av=repAB,ap=0))

      repL=list(av=1,ap=1)
      if(!("cAGE" %in% covariate[[t]][[meth]]$L)){
        repL$av = 0
      }
      if(!("SEX" %in% covariate[[t]][[meth]]$L)){
        repL$ap = 0
      }
      rep <- rbind(rep,data.frame(TypeOfSim=t,Method=meth,par="L",av=repL$av,ap=repL$ap))

      repS = 2
      if(!("SEX" %in% covariate[[t]][[meth]]$S) && !("RACE" %in% covariate[[t]][[meth]]$S)){
        repS = 0
      }else if(!("SEX" %in% covariate[[t]][[meth]]$S) || !("RACE" %in% covariate[[t]][[meth]]$S)){
        repS = 1
      }
      rep <- rbind(rep,data.frame(TypeOfSim=t,Method=meth,par="S",av=0,ap=repS))
    }
  }

  ## Selection In Model

  valueDisplay = data.frame()
  for(t in c(covType,corcovType)){
    for(meth in buildMethod){
      value=c(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                    & resultCovariateParCov$Method==meth
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="cAGE","ProportionSelected"]*100,
              resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                    & resultCovariateParCov$Method==meth
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="RACE","ProportionSelected"]*100,
              resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                    & resultCovariateParCov$Method==meth
                                    & resultCovariateParCov$TypeOfSim==t
                                    & resultCovariateParCov$Covariate=="SEX","ProportionSelected"]*100,
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                               & resultCovariateParCov$Method==meth
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="cAGE","ProportionSelected"]*100),digits=2),
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                               & resultCovariateParCov$Method==meth
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="RACE","ProportionSelected"]*100),digits=2),
              round(mean(resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                               & resultCovariateParCov$Method==meth
                                               & resultCovariateParCov$TypeOfSim==t
                                               & resultCovariateParCov$Covariate!="SEX","ProportionSelected"]*100),digits=2))
      valuemax = c(max(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[S]"
                                             & resultCovariateParCov$Method==meth
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="cAGE","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm=TRUE),
                   max(resultCovariateParCov[resultCovariateParCov$Parameter=="varphi[L]"
                                             & resultCovariateParCov$Method==meth
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="RACE","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm=TRUE),
                   max(resultCovariateParCov[resultCovariateParCov$Parameter=="delta[Ab]"
                                             & resultCovariateParCov$Method==meth
                                             & resultCovariateParCov$TypeOfSim==t
                                             & resultCovariateParCov$Covariate!="SEX","ProportionSelected"][(lim[[stringr::str_remove(t,"2")]]-5):lim[[stringr::str_remove(t,"2")]]]*100,na.rm=TRUE))

      valueDisplay <- rbind(valueDisplay,data.frame(cov=c(rep("cov",3),rep("FP",3)),
                                                    Method = meth,
                                                    value=value,
                                                    type=t,
                                                    coor=c(2,3,4,valuemax),
                                                    Parameter=rep(c("varphi[S]","varphi[L]","delta[Ab]"),2),
                                                    distext=paste0(c("cAGE :","RACE :","SEX :",rep("",3)),value,"%")))
    }
  }


  valueModel = cbind(resultModelCov,Parameter=c("delta[Ab]"))


  cmd = list()
  cmdA=list()
  for(t in c(covType,corcovType)){
    aux=list()
    auxA=list()
    for(meth in buildMethod){
      cmdAB = paste0('c(rep(gr,rep[rep$par=="AB" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","av"]),fill.vec[3],
                rep(gr,length(covariate[["',t,'"]][["',meth,'"]]$AB)- 1 - sum(rep[rep$par=="AB" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')
      cmdL = paste0('c(rep(gr,rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","av"]),fill.vec[2],rep(gr,rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","ap"]),
           rep(gr,length(covariate[["',t,'"]][["',meth,'"]]$L) - 1 - sum(rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')

      cmdS = paste0('c(fill.vec[1],rep(gr,rep[rep$par=="S" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","ap"]),
           rep(gr,length(covariate[["',t,'"]][["',meth,'"]]$S)-1 - sum(rep[rep$par=="S" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')

      cmdAAB = paste0('c(rep(0.2,rep[rep$par=="AB" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","av"]),0.5,
                rep(0.2,length(covariate[["',t,'"]][["',meth,'"]]$AB)- 1 - sum(rep[rep$par=="AB" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')
      cmdAL = paste0('c(rep(0.2,rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","av"]),0.5,rep(0.2,rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","ap"]),
           rep(0.2,length(covariate[["',t,'"]][["',meth,'"]]$L) - 1 - sum(rep[rep$par=="L" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')

      cmdAS = paste0('c(0.5,rep(0.2,rep[rep$par=="S" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'","ap"]),
           rep(0.2,length(covariate[["',t,'"]][["',meth,'"]]$S)-1 - sum(rep[rep$par=="S" & rep$TypeOfSim=="',t,'" & rep$Method=="',meth,'",c("av","ap")])))')


      aux <- append(aux,list(list(cmdAB=cmdAB,cmdL=cmdL,cmdS=cmdS)))
      names(aux)[length(aux)] = meth
      auxA <- append(auxA,list(list(cmdAAB=cmdAAB,cmdAL=cmdAL,cmdAS=cmdAS)))
      names(auxA)[length(auxA)] = meth
    }
    cmd <- append(cmd,list(aux))
    names(cmd)[length(cmd)] = t

    cmdA <- append(cmdA,list(auxA))
    names(cmdA)[length(cmdA)] = t
  }
  library(ggh4x)

  cmdCov = paste0(unlist(cmd$cov),collapse=",")
  eval(parse(text=paste0("fillScale=c(",cmdCov,")")))
  cmdACov = paste0(unlist(cmdA$cov),collapse=",")
  eval(parse(text=paste0("alphaScale=c(",cmdACov,")")))



  cov = ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==covType,],aes(x=factor(Covariate,levels=orderList$cov)))+
    geom_bar(position=position_dodge(preserve = "single"),alpha=alphaScale,
             color=fillScale,
             fill = fillScale)+
    facet_nested(factor(Method,levels=c("StepAIC","Lasso","ClustOfVar","Lasso with\nclustering step."))+Parameter~.,labeller=labeller(Parameter=label_parsed))+
    xlab("Covariate")+
    ylab("Count")+
    ggtitle("With uncorrelated covariates, ")+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    ylim(c(0,100))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=70,hjust=0),color=rep(c("#468b97","#ef6262", "#74C385"),length(buildMethod)),size=5.3)+
    geom_text(data=valueDisplay[valueDisplay$type==covType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$cov,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5.3)+
    geom_text(data=valueModel[valueModel$TypeOfSim==covType,],x=lim$cov,y=100,hjust=1,size=5,vjust=1,mapping=aes(label=text))+
    coord_cartesian(ylim=c(0,100),clip="off")


  cmdCov = paste0(unlist(cmd$corcov),collapse=",")
  eval(parse(text=paste0("fillScale=c(",cmdCov,")")))
  cmdACov = paste0(unlist(cmdA$corcov),collapse=",")
  eval(parse(text=paste0("alphaScale=c(",cmdACov,")")))


  corcov =    ggplot(CovariateModelSelectionCov[CovariateModelSelectionCov$TypeOfSim==corcovType,],aes(x=factor(Covariate,levels=orderList$corcov)))+
    geom_bar(position=position_dodge(preserve = "single"),alpha=alphaScale,
             color=fillScale,
             fill = fillScale)+
    facet_nested(factor(Method,levels=c("StepAIC","Lasso","ClustOfVar","Lasso with\nclustering step."))+Parameter~.,labeller=labeller(Parameter=label_parsed))+
    xlab("Covariate")+
    ylab("Count")+
    ggtitle("With correlated covariates, ")+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    theme(axis.text.y = element_text(size = 8))+
    ylim(c(0,100))+
    theme(axis.title = element_text(size=14))+
    theme(strip.text = element_text(size = 16))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.text = element_text(size=14))+
    theme(legend.title = element_text(size=12))+
    theme(plot.title = element_text(size=16,color="#ee6c4d"))+
    theme(plot.subtitle = element_text(size=12))+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="cov",],mapping=aes(label=distext,x=coor,y=70,hjust=0),color=rep(c("#468b97","#ef6262", "#74C385"),length(buildMethod)),size=5.3)+
    geom_text(data=valueDisplay[valueDisplay$type==corcovType & valueDisplay$cov=="FP",],mapping=aes(label=distext,x=lim$cov,y=coor+10,hjust=1,vjust=0),color="#9E9FA5",fontface = 'italic',size=5.3)+
    geom_text(data=valueModel[valueModel$TypeOfSim==corcovType,],x=lim$corcov,y=100,hjust=1,size=5,vjust=1,mapping=aes(label=text))+
    coord_cartesian(ylim=c(0,100),clip="off")

  annotate_figure(
    annotate_figure(ggarrange(cov, corcov, nrow = 2),
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".png"),
           height = 3000 + 1000*(length(buildMethod)-1), width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }


  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),".jpeg"),
           height = 3000+ 1000*(length(buildMethod)-1), width = 5500, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(cov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.png"),
           height = 3000 + 1000*(length(buildMethod)-1), width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }


  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Cov.jpeg"),
           height = 3000+ 1000*(length(buildMethod)-1), width = 5500, units = "px",device=grDevices::jpeg)
  }

  annotate_figure(
    annotate_figure(corcov,
                    top=text_grob(subtitle),
    ),
    top=text_grob("Covariate Selection Frequency",
                  face="bold",size=20,color="#862B0D")
  )


  if(PNG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.png"),
           height = 3000 + 1000*(length(buildMethod)-1), width = 5500, units = "px", bg='transparent',device=grDevices::png)
  }


  if(JPEG){
    ggsave(paste0(Folder,"/NumberSelectionParameter",covariateSize,"_",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),"_Corcov.jpeg"),
           height = 3000+ 1000*(length(buildMethod)-1), width = 5500, units = "px",device=grDevices::jpeg)
  }
}

######################################################################################################################################################

