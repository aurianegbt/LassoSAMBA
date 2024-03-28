tableStats <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Data.frame to use
  resultModelCov <- resultModel[resultModel$ProjectNumber == paste(covariateSize,"covariates") & resultModel$Method==buildMethod, ]
  
  resultModelParCov <- resultModelPar[resultModelPar$ProjectNumber == paste(covariateSize,"covariates") & resultModelPar$Method %in% buildMethod, ]

  errorStatsCov <- errorStats[errorStats$ProjectNumber == paste(covariateSize,"covariates") & errorStats$Method==buildMethod,]
  errorStatsParCov <- errorStatsPar[errorStatsPar$ProjectNumber == paste(covariateSize,"covariates") & errorStatsPar$Method==buildMethod,]
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
    mutate(FNR = (FN/sapply(FN+TN,FUN=function(x){max(x,1)})),.after="TN") %>%
    mutate(F1_score = TP/(TP+1/2*(FN+FP)),.after="FNR")
  

  df = data.frame(TypeOfSim=c("cov","corcov"),
                  FDR_mean = c(mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"]),
                               mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"])),
                  FDR_median = c(median(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"]),
                                 median(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"])),
                  FDR_sd = c(sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"]),
                             sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"])),
                  FDR_q975 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],0.975),
                               quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],0.975)),
                  FDR_q25 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FDR"],0.025),
                              quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FDR"],0.025)),
                  FNR_mean = c(mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"]),
                               mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"])),
                  FNR_median = c(median(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"]),
                                 median(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"])),
                  FNR_sd = c(sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"]),
                             sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"])),
                  FNR_q975 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],0.975),
                               quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],0.975)),
                  FNR_q25 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FNR"],0.025),
                              quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FNR"],0.025)),
                  FN_mean = c(mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FN"]),
                              mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FN"])),
                  FN_median = c(median(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FN"]),
                                median(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FN"])),
                  FN_sd = c(sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FN"]),
                            sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FN"])),
                  FN_q975 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FN"],0.975),
                              quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FN"],0.975)),
                  FN_q25 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","FN"],0.025),
                             quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","FN"],0.025)),
                  F1_score_mean = c(mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","F1_score"]),
                              mean(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","F1_score"])),
                  F1_score_median = c(median(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","F1_score"]),
                                median(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","F1_score"])),
                  F1_score_sd = c(sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","F1_score"]),
                            sd(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","F1_score"])),
                  F1_score_q975 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","F1_score"],0.975),
                              quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","F1_score"],0.975)),
                  F1_score_q25 = c(quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="cov","F1_score"],0.025),
                             quantile(errorStatsParCov[errorStatsParCov$TypeOfSim=="corcov","F1_score"],0.025)))
  # Function
  percent <- function(x, digits = 1, format = "f", ...) {
    paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
  }

  CB <- function(df,char){
    lb = df[,paste0(char,"_q25")]
    ub = df[,paste0(char,"_q975")]
    return(paste0("[",sapply(lb,FUN=function(x){percent(max(0,x))}),
                  ";",sapply(ub,FUN=function(x){percent(min(1,x))}),"]"))}

  # Data processing
  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"),FN_CB=CB(df,"FN"),F1_score_CB=CB(df,"F1_score"))

  tableCov = data.frame(c("With uncorrelated covariates :", "False Discovery Rate","False Negative Rate","F1 score",
                          paste0("\t • Final Final model without any False Negatives : ",
                                 percent(resultModelParCov[resultModelParCov$TypeOfSim=="cov","NoFNModel"],digits=0)),
                          paste0("\t • Final model is the true one : ",
                                 percent(resultModelParCov[resultModelParCov$TypeOfSim=="cov","TrueModel"],digits=0))),

                        median=c("",percent((as.numeric(format(c(df[df$TypeOfSim=="cov","FDR_median"],
                                                                 df[df$TypeOfSim=="cov","FNR_median"],
                                                                 df[df$TypeOfSim=="cov","F1_score_median"]), scientific=TRUE)))),"",""),
                        CB=c("",df[df$TypeOfSim=="cov","FDR_CB"],
                             df[df$TypeOfSim=="cov","FNR_CB"],
                             df[df$TypeOfSim=="cov","F1_score_CB"],"",""))

  colnames(tableCov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  tableCorcov =data.frame(c("With correlated covariates :","False Discovery Rate","False Negative Rate","F1 score",
                            paste0("\t • Final Final model without any False Negatives : ",
                                   percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"],digits=0)),
                            paste0("\t • Final model is the true one : ",
                                   percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"],digits=0))),

                          median=c("",percent(as.numeric(format(c(df[df$TypeOfSim=="corcov","FDR_median"],
                                                                  df[df$TypeOfSim=="corcov","FNR_median"],
                                                                  df[df$TypeOfSim=="corcov","F1_score_median"]), scientific=TRUE))),"",""),
                          CB=c("",df[df$TypeOfSim=="corcov","FDR_CB"],df[df$TypeOfSim=="corcov","FNR_CB"],df[df$TypeOfSim=="corcov","F1_score_CB"],"",""))

  colnames(tableCorcov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  table <-  tibble::as_tibble(rbind(tableCov,tableCorcov))


  # Table for stats info

  ft <- flextable(table) %>%
    merge_at(i=1,j=1:3) %>%
    merge_at(i=5,j=1:3) %>%
    merge_at(i=6,j=1:3) %>%
    merge_at(i=7,j=1:3) %>%
    merge_at(i=11,j=1:3) %>%
    merge_at(i=12,j=1:3) %>%
    color(i=c(1,7),color="indianred") %>%
    bold(i=c(1,7)) %>%
    set_table_properties(layout="autofit",width=1) %>%
    bg(i=c(5:6,11:12),j=1:3,bg="#e4e6eb",part="body") %>%
    hline(i = 6, part = "body", border = fp_border_default(color = "grey1", width = 1) ) %>%
    hline(i = 1, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    hline(i = 7, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    vline(j=c(2,1),i=c(2,3,4,8,9,10),part="body", border = fp_border_default(color = "grey", width = 1)) %>%
    add_header_lines("Error Rate Comparison Table") %>%
    bold(part="header") %>%
    fontsize(size=22,part="header",i=1)%>%
    fontsize(size=20,part="header",i=2)%>%
    fontsize(size=18,part="body") %>%
    fontsize(size=20,i=c(1,7),part="body") %>%
    align(align = "center", part = "header",i=1) %>%
    # fontsize(size=16,i=c(6,7,11,12),part="body") %>%
    add_footer_lines("Covariates presence in final model with parameters link")%>%
    add_footer_lines(subtitle) %>%
    fontsize(size=16,part="footer") %>%
    italic(i=1,par="footer")%>%
    fontsize(i=1,size=12,part="footer") %>%
    align( i =1, align="right",part="footer") %>%
    align( i =2, align="left",part="footer")
    

  # Save plot
  save_as_html(ft, path = paste0(Folder,"/ErrorTable.html"),expand=10)
  webshot(url=paste0(Folder,"/ErrorTable.html"),file=paste0(Folder,"/ErrorTable.png"),quiet=TRUE)
  unlink(paste0(Folder,"/ErrorTable.html"))
}
