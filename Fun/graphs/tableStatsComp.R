tableStatsComp <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Data.frame to use
  resultModelCov <- resultModel[resultModel$ProjectNumber == paste(covariateSize,"covariates") & resultModel$Method %in% buildMethod, ]
  
  
  resultModelParCov <- resultModelPar[resultModelPar$ProjectNumber == paste(covariateSize,"covariates") & resultModelPar$Method %in% buildMethod, ]

  errorStatsCov <- errorStats[errorStats$ProjectNumber == paste(covariateSize,"covariates") & errorStats$Method %in% buildMethod,]
  errorStatsParCov <- errorStatsPar[errorStatsPar$ProjectNumber == paste(covariateSize,"covariates") & errorStatsPar$Method %in% buildMethod,]
  errorStatsParCov <-  suppressMessages(errorStatsParCov %>%
    group_by(Model,ProjectNumber,TypeOfSim,Method) %>%
    summarise(
      FP = sum(FP),
      TP = sum(TP), 
      FN = sum(FN),
      TN = sum(TN)
    )) %>%
    as.data.frame()

  errorStatsParCov <- errorStatsParCov  %>% 
    mutate(FDR = (FP/sapply(FP+TP,FUN=function(x){max(x,1)})),.after = "TP") %>%
    mutate(FNR = (FN/sapply(FP+TP,FUN=function(x){max(x,1)})),.after="TN") 

  df = data.frame(TypeOfSim=rep(c("cov","corcov"),each=length(buildMethod)),
                  Method = rep(buildMethod,2),
                  FDR_mean = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN= mean)[buildMethod],
                               sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = mean)[buildMethod]),

                  FDR_median = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],
                                              errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN = median)[buildMethod],
                                 sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],
                                              errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = median)[buildMethod]),

                  FDR_sd = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],
                                          errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN = sd)[buildMethod],
                             sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],
                                          errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = sd)[buildMethod]),

                  FDR_q975 = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),
                                      FUN = function(x){return(quantile(x,0.975))})[paste0(buildMethod,".97.5%")],
                               sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),
                                      FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")]),

                  FDR_q25 = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],
                                           errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")],
                              sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],
                                           errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")]),

                  FNR_mean = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN = mean)[buildMethod],
                               sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],
                                            errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = mean)[buildMethod]),

                  FNR_median = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],
                                              errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN = median)[buildMethod],
                                 sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],
                                              errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = median)[buildMethod]),

                  FNR_sd = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],
                                          errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),FUN = sd)[buildMethod],
                             sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],
                                          errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),FUN = sd)[buildMethod]),

                  FNR_q975 =  c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],
                                             errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),
                                       FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")],
                                sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],
                                             errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),
                                       FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")]),

                  FNR_q25 = c(sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],
                                           errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","Method"]),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")],
                              sapply(split(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],
                                           errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","Method"]),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")]))

  # Function
  percent <- function(x, digits = 1, format = "f", ...) {      # Create user-defined function
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
  }

  CB <- function(df,char){
    lb = df[,paste0(char,"_q25")]
    ub = df[,paste0(char,"_q975")]
    return(paste0("[",sapply(lb,FUN=function(x){percent(max(0,x))}),
                  ";",sapply(ub,FUN=function(x){percent(min(1,x))}),"]"))}

  # Data Processing
  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"))
  methname=c(reg="stepAIC",lasso="Lasso\nwhithout s.s.",elastinet="Elastic net\nwhithout s.s.",lassoSS="Lasso",elasticnetSS="Elastic Net",rlasso="Lasso with\ns.s. on replicates",relasticnet="Elastic Net with\ns.s. on replicates",lassoCrit = "Lasso with\nmultiple thresholds and no s.s.",elasticnetCrit="Elastic Net with\nmultiple thresholds and no s.s.",lassoSSCrit="Lasso with\nmultiple thresholds",elasticnetSSCrit="Elastic Net with\nmultiple thresholds",rlassoCrit="Lasso with mult.\nthresholds and s.s. on rep.",relasticnetCrit="Elastic Net with mult.\nthresholds and s.s. on rep.",setNames(paste0("penalized stepAIC\npen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),regnoCov0="stepAIC\nwhithout stat. test",lassoSSCov0="Lasso\nwhithout stat. test",elasticnetSSnoCov0="Elastic Net\nwhithout stat. test")

  tableCov = data.frame(Rate = c("With uncorrelated covariates :",
                                 paste0("False Discovery Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0("False Negative Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0(" • Final Final model without any False Negatives :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelParCov[resultModelParCov$TypeOfSim=="cov","NoFNModel"],resultModelParCov[resultModelParCov$TypeOfSim=="cov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n")),
                                 paste0("  • Final model is the true one :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelParCov[resultModelParCov$TypeOfSim=="cov","TrueModel"],resultModelParCov[resultModelParCov$TypeOfSim=="cov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n"))),

                        Median = c("",paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="cov",c("FDR_median")],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],percent),collapse="\n")),
                                   paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="cov",c("FNR_median")],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],percent),collapse="\n")),"",""),

                        CB = c("",paste0("\n",paste0(split(df[df$TypeOfSim=="cov","FDR_CB"],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],collapse="\n")),
                               paste0("\n",paste0(split(df[df$TypeOfSim=="cov","FNR_CB"],df[df$TypeOfSim=="cov",c("Method")])[buildMethod],collapse="\n")),"",""))

  colnames(tableCov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  tableCorcov = data.frame(Rate = c("With correlated covariates :",
                                    paste0("False Discovery Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                    paste0("False Negative Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                    paste0("• Final Final model without any False Negatives :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelParCov[resultModelParCov$TypeOfSim=="corcov","NoFNModel"],resultModelParCov[resultModelParCov$TypeOfSim=="corcov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n")),
                                    paste0("• Final model is the true one :\n ",paste0(paste0("\t\t - ",methname[buildMethod]," : ",sapply(split(resultModelParCov[resultModelParCov$TypeOfSim=="corcov","TrueModel"],resultModelParCov[resultModelParCov$TypeOfSim=="corcov","Method"])[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n"))),

                           Median = c("",paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="corcov",c("FDR_median")],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],percent),collapse="\n")),
                                      paste0("\n",paste0(sapply(split(df[df$TypeOfSim=="corcov",c("FNR_median")],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],percent),collapse="\n")),"",""),

                           CB = c("",paste0("\n",paste0(split(df[df$TypeOfSim=="corcov","FDR_CB"],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],collapse="\n")),
                                  paste0("\n",paste0(split(df[df$TypeOfSim=="corcov","FNR_CB"],df[df$TypeOfSim=="corcov",c("Method")])[buildMethod],collapse="\n")),"",""))

  colnames(tableCorcov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")


  table <-  tibble::as_tibble(rbind(tableCov,tableCorcov))

  # Table with stats info
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
    fontsize(size=16,i=c(4,5,9,10),part="body")%>%
    add_footer_lines("Covariates presence in final model with parameters link")%>%
    add_footer_lines(subtitle) %>%
    fontsize(size=16,part="footer") %>%
    italic(i=1,par="footer")%>%
    fontsize(i=1,size=12,part="footer") %>%
    align( i =1, align="right",part="footer") %>%
    align( i =2, align="left",part="footer")

  # Save plot
  save_as_html(ft, path = paste0(Folder,"/ErrorTable",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),covariateSize,".html"),expand=10)
  
    webshot(paste0(Folder,"/ErrorTable",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),covariateSize,".html"), paste0(Folder,"/ErrorTable",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),covariateSize,".png"),quiet=TRUE)
    
  unlink(paste0(Folder,"/ErrorTable",paste0(sapply(buildMethod,function(x){toupper(stringr::str_sub(x,end=2))}),collapse="-"),covariateSize,".html"))
}
