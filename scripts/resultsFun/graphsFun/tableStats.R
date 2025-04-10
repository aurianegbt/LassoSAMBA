tableStats <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))

  # Data.frame to use
  resultModelParCov <- resultModelPar[resultModelPar$Method==buildMethod, ]
  errorStatsParCov <- errorStatsPar[errorStatsPar$Method==buildMethod,]
  
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
    mutate(FPR = (FP/sapply(FP+TN,FUN=function(x){max(x,1)})),.after = "TP") %>%
    mutate(FDR = (FP/sapply(FP+TP,FUN=function(x){max(x,1)})),.after = "TP") %>%
    mutate(FNR = (FN/sapply(FN+TP,FUN=function(x){max(x,1)})),.after="TN") %>%
    mutate(FOR = (FN/sapply(FN+TN,FUN=function(x){max(x,1)})),.after="TN") %>%
    mutate(F1_score = TP/(TP+1/2*(FN+FP)),.after="FNR")
  

  df = data.frame(FDR_mean = mean(errorStatsParCov[,"FDR"]),
         FDR_median = median(errorStatsParCov[,"FDR"]),
         FDR_sd = sd(errorStatsParCov[,"FDR"]),
         FDR_q975 = quantile(errorStatsParCov[,"FDR"],0.975),
         FDR_q25 = quantile(errorStatsParCov[,"FDR"],0.025),
         
         FNR_mean = mean(errorStatsParCov[,"FNR"]),
         FNR_median = median(errorStatsParCov[,"FNR"]),
         FNR_sd = sd(errorStatsParCov[,"FNR"]),
         FNR_q975 = quantile(errorStatsParCov[,"FNR"],0.975),
         FNR_q25 = quantile(errorStatsParCov[,"FNR"],0.025),
         
         
         FOR_mean = mean(errorStatsParCov[,"FOR"]),
         FOR_median = median(errorStatsParCov[,"FOR"]),
         FOR_sd = sd(errorStatsParCov[,"FOR"]),
         FOR_q975 = quantile(errorStatsParCov[,"FOR"],0.975),
         FOR_q25 = quantile(errorStatsParCov[,"FOR"],0.025),
         
         
         FPR_mean = mean(errorStatsParCov[,"FPR"]),
         FPR_median = median(errorStatsParCov[,"FPR"]),
         FPR_sd = sd(errorStatsParCov[,"FPR"]),
         FPR_q975 = quantile(errorStatsParCov[,"FPR"],0.975),
         FPR_q25 = quantile(errorStatsParCov[,"FPR"],0.025),
         
         FN_mean = mean(errorStatsParCov[,"FN"]),
         FN_median = median(errorStatsParCov[,"FN"]),
         FN_sd = sd(errorStatsParCov[,"FN"]),
         FN_q975 = quantile(errorStatsParCov[,"FN"],0.975),
         FN_q25 = quantile(errorStatsParCov[,"FN"],0.025),
         
         F1_score_mean = mean(errorStatsParCov[,"F1_score"]),
         F1_score_median = median(errorStatsParCov[,"F1_score"]),
         F1_score_sd = sd(errorStatsParCov[,"F1_score"]),
         F1_score_q975 = quantile(errorStatsParCov[,"F1_score"],0.975),
         F1_score_q25 = quantile(errorStatsParCov[,"F1_score"],0.025))
  
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

  table = data.frame(c("False Discovery Rate","False Negative Rate","F1 score",
                          paste0("\t • Final Final model without any False Negatives : ",
                                 percent(resultModelParCov[,"NoFNModel"],digits=0)),
                          paste0("\t • Final model is the true one : ",
                                 percent(resultModelParCov[,"TrueModel"],digits=0))),

                        median=c(percent((as.numeric(format(c(df[,"FDR_median"],
                                                                 df[,"FNR_median"],
                                                                 df[,"F1_score_median"]), scientific=TRUE)))),"",""),
                        CB=c(df[,"FDR_CB"],
                             df[,"FNR_CB"],
                             df[,"F1_score_CB"],"",""))

  colnames(table) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")
  
  table <-  tibble::as_tibble(table)


  # Table for stats info
  
  # writeLines(paste0(table$Median,"  ",table$`Confidence Interval\n(quantiles 95%)`),con=paste0(Folder,"/ErrorTable.txt"))

  ft <- flextable(table) %>%
    merge_at(i=4,j=1:3) %>%
    merge_at(i=5,j=1:3) %>%
    set_table_properties(layout="autofit",width=1) %>%
    bg(i=4:5,j=1:3,bg="#e4e6eb",part="body") %>%
    hline(i = 5, part = "body", border = fp_border_default(color = "grey1", width = 1) ) %>%
    hline(i = 1, part = "body", border = fp_border_default(color = "grey", width = 0.7) ) %>%
    vline(j=1:2,i=1:3,part="body", border = fp_border_default(color = "grey", width = 1)) %>%
    add_header_lines("Error Rate Comparison Table") %>%
    bold(part="header") %>%
    color(i=1,part="header",color = "indianred4") %>%
    color(i=2,part="header",color = "indianred") %>%
    fontsize(size=24,part="header",i=1)%>%
    fontsize(size=20,part="header",i=2)%>%
    fontsize(size=18,part="body") %>%
    align(align = "center", part = "header",i=1) %>%
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
