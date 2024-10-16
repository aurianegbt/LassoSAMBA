tableStatsComp <- function(Folder,subtitle,project,buildMethod,JPEG,PNG){
  
  # Load data
  load(paste0("outputs/finalResults/BuildResults_",project,".RData"))
  source(paste0("data/simulationFiles/Files",project,"/H1.all.R"))
  

  # Data.frame to use
  resultModelCov <- resultModel[resultModel$Method %in% buildMethod, ]
  
  
  resultModelParCov <- resultModelPar[resultModelPar$Method %in% buildMethod, ]

  errorStatsCov <- errorStats[errorStats$Method %in% buildMethod,]
  errorStatsParCov <- errorStatsPar[errorStatsPar$Method %in% buildMethod,]
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
    mutate(FDR = (FP/sapply(FP+TP,FUN=function(x){max(x,1)})),.after = "TP") %>%
    mutate(FNR = (FN/sapply(FN+TN,FUN=function(x){max(x,1)})),.after="TN") %>%
    mutate(F1_score = TP/(TP+1/2*(FN+FP)),.after="FNR")

  df = data.frame(Method = buildMethod,
                  FDR_mean = sapply(split(errorStatsParCov$FDR,
                                            errorStatsParCov$Method),FUN= mean)[buildMethod],

                  FDR_median = sapply(split(errorStatsParCov$FDR,
                                              errorStatsParCov$Method),FUN = median)[buildMethod],

                  FDR_sd = sapply(split(errorStatsParCov$FDR,
                                          errorStatsParCov$Method),FUN = sd)[buildMethod],

                  FDR_q975 = sapply(split(errorStatsParCov$FDR,
                                            errorStatsParCov$Method),
                                      FUN = function(x){return(quantile(x,0.975))})[paste0(buildMethod,".97.5%")],

                  FDR_q25 = sapply(split(errorStatsParCov$FDR,
                                           errorStatsParCov$Method),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")],

                  FNR_mean = sapply(split(errorStatsParCov$FNR,
                                            errorStatsParCov$Method),FUN = mean)[buildMethod],

                  FNR_median = sapply(split(errorStatsParCov$FNR,
                                              errorStatsParCov$Method),FUN = median)[buildMethod],
                                 
                  FNR_sd = sapply(split(errorStatsParCov$FNR,
                                          errorStatsParCov$Method),FUN = sd)[buildMethod],
                             
                  FNR_q975 =  sapply(split(errorStatsParCov$FNR,
                                             errorStatsParCov$Method),
                                       FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")],

                  FNR_q25 = sapply(split(errorStatsParCov$FNR,
                                           errorStatsParCov$Method),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")],
                  
                  F1_score_mean = sapply(split(errorStatsParCov$F1_score,
                                            errorStatsParCov$Method),FUN = mean)[buildMethod],
                  
                  F1_score_median = sapply(split(errorStatsParCov$F1_score,
                                              errorStatsParCov$Method),FUN = median)[buildMethod],
                  
                  F1_score_sd = sapply(split(errorStatsParCov$F1_score,
                                          errorStatsParCov$Method),FUN = sd)[buildMethod],
                  
                  F1_score_q975 =  sapply(split(errorStatsParCov$F1_score,
                                             errorStatsParCov$Method),
                                       FUN = function(x){quantile(x,0.975)})[paste0(buildMethod,".97.5%")],
                  
                  F1_score_q25 = sapply(split(errorStatsParCov$F1_score,
                                           errorStatsParCov$Method),
                                     FUN = function(x){quantile(x,0.025)})[paste0(buildMethod,".2.5%")])

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
  df <- cbind(df, FDR_CB = CB(df,"FDR"),FNR_CB = CB(df,"FNR"),F1_score_CB=CB(df,"F1_score"))
  methname=c(reg="stepAIC with stat. test",
             lasso="Lasso without s.s.",
             elastinet="Elastic net without s.s.",
             lassoSS="Lasso with stat. test",
             elasticnetSS="Elastic Net with stat. test",
             rlasso="Lasso with s.s. on replicates",
             relasticnet="Elastic Net with s.s. on replicates",
             lassoCrit = "Lasso with multiple thresholds and no s.s.",
             elasticnetCrit="Elastic Net with multiple thresholds and no s.s.",
             lassoSSCrit="Lasso with multiple thresholds",
             elasticnetSSCrit="Elastic Net with multiple thresholds",
             rlassoCrit="Lasso with mult. thresholds and s.s. on rep.",
             relasticnetCrit="Elastic Net with mult. thresholds and s.s. on rep.",
             setNames(paste0("penalized stepAIC pen=",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regPEN")),buildMethod[stringr::str_detect(buildMethod,"regPEN")]),
             regnoCov0="stepAIC",
             lassoSSnoCov0="Lasso",
             lassoSSCritnoCov0="Lasso",
             SAEMVS="SAEMVS",
             setNames(paste0("Lasso ",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))],"sharpnoCov0"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0"))]),
             setNames(paste0("Lasso : E[FDR]<",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))],"sharpnoCov0FDP"),"%"),buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharpnoCov0FDP"))]),
             setNames(paste0("Lasso ",stringr::str_remove_all(buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))],"sharp"),"% higher score"),buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove(buildMethod,"sharp"))]),
             elasticnetSSnoCov0="Elastic Net",
             sharpnoCov0="Lasso calibrated using sharp",
             sharp="Lasso calibrated using sharp with stat. test")

  table = data.frame(Rate = c(paste0("False Discovery Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0("False Negative Rate :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0("F1 score :\n",paste0(paste0("\t\t - ",methname[buildMethod]),collapse="\n")),
                                 paste0(" • Final Final model without any False Negatives :\n ",paste0(paste0("\t\t - ",stringr::str_replace_all(methname[buildMethod],"\n"," ")," : ",sapply(split(resultModelParCov[,"NoFNModel"],resultModelParCov$Method)[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n")),
                                 paste0("  • Final model is the true one :\n ",paste0(paste0("\t\t - ",stringr::str_replace_all(methname[buildMethod],"\n"," ")," : ",sapply(split(resultModelParCov$TrueModel,resultModelParCov$Method)[buildMethod],FUN=function(x){percent(x,digits=0)})[buildMethod]),collapse="\n"))),

                        Median = c(paste0("\n",paste0(sapply(split(df$FDR_median,df$Method)[buildMethod],percent),collapse="\n")),
                                   paste0("\n",paste0(sapply(split(df$FNR_median,df$Method)[buildMethod],percent),collapse="\n")),
                                   paste0("\n",paste0(sapply(split(df$F1_score_median,df$Method)[buildMethod],percent),collapse="\n")),"",""),

                        CB = c(paste0("\n",paste0(split(df$FDR_CB,df$Method)[buildMethod],collapse="\n")),
                               paste0("\n",paste0(split(df$FNR_CB,df$Method)[buildMethod],collapse="\n")),
                               paste0("\n",paste0(split(df$F1_score_CB,df$Method)[buildMethod],collapse="\n")),"",""))

  colnames(table) <- c("Rate","Median","Confidence Interval\n(quantiles 95%)")
  
  table <-  tibble::as_tibble(table)

  writeLines(c(paste0(unlist(strsplit(table$Median,"\n")),"  ",unlist(strsplit(table$`Confidence Interval\n(quantiles 95%)`,"\n"))),"  ",unname(sapply(split(resultModelParCov$NoFNModel,resultModelParCov$Method)[buildMethod],FUN=function(x){percent(x,digits=0)})),"  ",unname(sapply(split(resultModelParCov$TrueModel,resultModelParCov$Method)[buildMethod],FUN=function(x){percent(x,digits=0)}))),con=paste0(Folder,"/ErrorTable.txt"))
  
  # Table with stats info
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
  
    webshot(paste0(Folder,"/ErrorTable.html"), paste0(Folder,"/ErrorTable.png"),quiet=TRUE)
    
  unlink(paste0(Folder,"/ErrorTable.html"))
}
