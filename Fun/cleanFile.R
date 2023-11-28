cleanFile <- function(project=c("Pasin","Warfarine"),
                      covariateSize=c(10,50,200,500),
                      covariateType=c("cov","corcov"),
                      buildMethod=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"),
                      completePrint=FALSE,
                      print=FALSE){
  for(proj in project){
    for(cov in covariateSize){
      for(type in covariateType){
        for(meth in buildMethod){
          
          resultsFolder = paste0("Results",proj,"/Results_",meth,"/",cov,type)
          
          deletedFiles <- numeric()
          for(i in 1:100){
            comp = paste0(resultsFolder,"/compResults_",i,".RData")
            build = paste0(resultsFolder,"/buildResults_",i,".RData")
            
            ## LE SEUL CAS OK : LES DEUX EXISTENT ET SONT NON VIDE 
            if(file.exists(comp) || file.exists(build)){
              errCOMP = !(fileTEST(comp))
              errBUILD = !(fileTEST(build))
              
              if(errCOMP || errBUILD){
                unlink(comp)
                unlink(build)
                deletedFiles <- append(deletedFiles,i)
              }
            }
          }
          if(print && length(deletedFiles)!=0){
            cat(paste0("▶ For project with ",cov))
            if(type=="corcov"){cat(" correlated covariates, ")}else{cat(" covariates,")}
            cat(" built with ")
            if(meth == "reg"){
              cat("stepAIC regression")
            }else if(meth=="lassoSS"){
              cat("lasso with stability selection")
            }else if(meth=="lassoSSCrit"){
              cat("lasso with multiple thresholds and stability selection")
            }else if(meth=="elasticnetSS"){
              cat("elastic net with stability selection")
            }else if(meth=="elasticnetSSCrit"){
              cat("elastic net with multiple thresholds and stability selection")
            }
            cat(" :")
            first = FALSE
            dFlen = length(deletedFiles)
            if(dFlen==1){
              cat(paste0("   • File n°",deletedFiles[1]," has been deleted.\n"))
            }else{
              cat(paste0("   • File n°",paste0(paste0(deletedFiles[-dFlen],collapse=", ")," and ",deletedFiles[dFlen])," have been deleted.\n"))
            }
          }
        }
      }
    } 
  }
}

  

fileTEST <- function(link){
    bool = tryCatch({
      load(link)
      bool <- TRUE
    }, error = function(e) {
      print(e)
      bool <- FALSE
    })
    return(bool)
}




