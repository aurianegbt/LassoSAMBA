#' Sumarize results obtain From Built model. 
#'
#' @param pathToRes 
#' @param Latex wether or not to display the latex code for the table.

resultsFB <- function(projectNameNumber,projectType,populationSize,criterion="BICc",PREVAc=FALSE,Latex=FALSE,meth="reg",opt="standard"){
  pathToRes = paste0("Results/Results_",opt,"/Results_",meth,"/finaleResults/",projectNameNumber,projectType,".RData")
  CovariateModelSelection <- data.frame() 
  
  load(pathToRes)
       
  n = length(Models)
  # allCovMentionned = character(0)
  for (j in 1:n){
    Model = Models[[j]]
    par = c("phi_S","phi_L","delta_AB")
    allCovMentionned = character(0)
    cov = names(Model$covariateModel[[1]])
    for(c in cov){if(!(c %in% allCovMentionned)){allCovMentionned <- c(allCovMentionned,c)}}
      
    if(length(names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]))!=0){auxS = data.frame(Covariate = names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]),
                        Parameter = "phi_S",Model=as.character(j))
    }else{ auxS = data.frame()}
    if(length(names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]))!=0){auxL = data.frame(Covariate = names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]),
                          Parameter = "phi_L",Model=as.character(j))
    }else{ auxL = data.frame()}
    if(length(names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]))!=0){ auxD =data.frame(Covariate = names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]),
                        Parameter = "delta_AB",Model=as.character(j))
     }else{auxD = data.frame()}
    CovariateModelSelection <- rbind(CovariateModelSelection,auxS,auxL,auxD)
  }
  
  resultCovariate1 <- data.frame()  # for AGE SEX RACE
  resultCovariate0 <- data.frame()  # other covariate 
  covLink = c(cAGE = "phi_S",RACE="phi_L",SEX="delta_AB")
  
  for(covariate in allCovMentionned){
    df = CovariateModelSelection[CovariateModelSelection$Covariate==covariate,]
    TotalNumberofModel = n
    TotalSelected = nrow(df)
    SelectedInDistinctModel = length(unique(df[,"Model"]))
      
    if(covariate %in% c("cAGE","RACE","SEX")){
      MissingInModel = TotalNumberofModel - SelectedInDistinctModel
      WellDefinedInModel = 0
      for(model in unique(df[,"Model"])){
        if(nrow(df[df$Model==model,])==1 && df[df$Model==model,"Parameter"]==covLink[covariate]){
          WellDefinedInModel = WellDefinedInModel + 1
        }
      } 
      resultCovariate1 <- rbind(resultCovariate1,data.frame(Covariate =covariate,TotalNumberofModel,TotalSelected,
                                                            SelectedInDistinctModel,MissingInModel,WellDefinedInModel))
    }else{
        resultCovariate0 <- rbind(resultCovariate0,data.frame(Covariate =covariate,TotalNumberofModel,TotalSelected,SelectedInDistinctModel))
    }
  }
  
  if(Latex){
    cat("\n \n ----------- LATEX CODE FOR RESULTS -----------\n")
    colnames(resultCovariate1)[-c(1,2)] <- c("Number of Selection","Selected in Model (%)","Missing in Model (%)","Correct Model (%)")
    tot = unique(resultCovariate1[,"TotalNumberofModel"])
    df = resultCovariate1[,-c(2,5)]
    df[,-c(1:2)] <- df[,-c(1:2)]/ tot
    print(xtable::xtable(df, digits  = 2),include.rownames=FALSE)
    cat("\n ----------------------------------------------\n")
    colnames(resultCovariate0)[c(3,4)] <- c("Number of Selection","Selected in Model (%)")
    tot = unique(resultCovariate0[,"TotalNumberofModel"])
    df = resultCovariate0[,-2]
    df[,3] <- df[,3]/ tot
    print(xtable::xtable(df, digits  = 2),include.rownames=FALSE)
    cat("\n ----------------------------------------------\n")
  }
  if(Res)
    return(list(CovariateModelSelection,resultCovariate1,resultCovariate0))  
}