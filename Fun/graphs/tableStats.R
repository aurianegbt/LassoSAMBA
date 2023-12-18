tableStats <- function(Folder,subtitle,project,covariateSize,buildMethod,JPEG,PNG){

  # Load data
  load(paste0("Save/BuildResults_",project,".RData"))
  source(paste0("Files/Files",project,"/H1.all.R"))

  # Data.frame to use
  resultModelCov <- resultModel[resultModel$ProjectNumber == paste(covariateSize,"covariates") & resultModel$Method==buildMethod, ]

  errorStatsCov <- errorStats[errorStats$ProjectNumber == paste(covariateSize,"covariates") & errorStats$Method==buildMethod,]

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
  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"),FN_CB=CB(df,"FN"))

  tableCov = data.frame(c("With uncorrelated covariates :", "False Discovery Rate","False Negative Rate",
                          paste0("\t • Final Final model without any False Negatives : ",
                                 percent(resultModelCov[resultModelCov$TypeOfSim=="cov","NoFNModel"],digits=0)),
                          paste0("\t • Final model is the true one : ",
                                 percent(resultModelCov[resultModelCov$TypeOfSim=="cov","TrueModel"],digits=0))),

                        median=c("",percent((as.numeric(format(c(df[df$TypeOfSim=="cov","FDR_median"],
                                                                 df[df$TypeOfSim=="cov","FNR_median"])), scientific=TRUE, digits=2))),"",""),
                        CB=c("",df[df$TypeOfSim=="cov","FDR_CB"],
                             df[df$TypeOfSim=="cov","FNR_CB"],"",""))

  colnames(tableCov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  tableCorcov =data.frame(c("With correlated covariates :","False Discovery Rate","False Negative Rate",
                            paste0("\t • Final Final model without any False Negatives : ",
                                   percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","NoFNModel"],digits=0)),
                            paste0("\t • Final model is the true one : ",
                                   percent(resultModelCov[resultModelCov$TypeOfSim=="corcov","TrueModel"],digits=0))),

                          median=c("",percent(as.numeric(format(c(df[df$TypeOfSim=="corcov","FDR_median"],
                                                                  df[df$TypeOfSim=="corcov","FNR_median"]), scientific=TRUE, digits=2))),"",""),
                          CB=c("",df[df$TypeOfSim=="corcov","FDR_CB"],df[df$TypeOfSim=="corcov","FNR_CB"],"",""))

  colnames(tableCorcov) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")

  table <-  tibble::as_tibble(rbind(tableCov,tableCorcov))


  # Table for stats info

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
    add_footer_lines("Covariates presence in final model with parameters link")%>%
    add_footer_lines(subtitle) %>%
    fontsize(size=16,part="footer") %>%
    italic(i=1,par="footer")%>%
    fontsize(i=1,size=12,part="footer") %>%
    align( i =1, align="right",part="footer") %>%
    align( i =2, align="left",part="footer")
    

  # Save plot
  save_as_html(ft, path = paste0(Folder,"/ErrorTable.html"),expand=10)
  if(PNG){
    webshot(url=paste0(Folder,"/ErrorTable.html"),file=paste0(Folder,"/ErrorTable.png"),quiet=TRUE)
  }
  if(JPEG){
    webshot(url=paste0(Folder,"/ErrorTable.html"),file=paste0(Folder,"/ErrorTable.jpeg"),quiet=TRUE)
  }
  unlink(paste0(Folder,"/ErrorTable.html"))
}
