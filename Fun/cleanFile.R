cleanFile <- function(covariateSize=c(100,200),
                      covariateType=c("cov","corcov"),
                      completePrint=FALSE,
                      buildMethod=c("reg","lassoSS","lassoSSCl"),
                      print=TRUE,
                      buildOptions=c("standard"),
                      Rsmlx=c("","2","Step")){
  
  for(Rv in Rsmlx){
    for(opt in buildOptions){
      for(cov in covariateSize){
        
        for(type in covariateType){
          for(meth in buildMethod){
            
            resultsFolder =paste0("Results/Results",Rv,"/Results_",meth,"/",cov,type)
            
            for(resFD in resultsFolder){
              deletedFiles <- numeric()
              for(i in 1:150){
                comp = paste0(resFD,"/compResults_",i,".RData")
                build = paste0(resFD,"/buildResults_",i,".RData")
                
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
                }else if(meth=="lasso"){
                  cat("lasso regression")
                }else{
                  cat("lasso with stability selection")
                }
                if(opt=="noCov0"){
                  cat(", whithout including covariates :")
                }else if(opt=="Mstar"){
                  cat(", starting from M* model :")
                }else{
                  cat(" :")
                }
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




