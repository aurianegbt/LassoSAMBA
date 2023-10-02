gatherResults <- function(covariateSize=c(100,200),
                          covariateType=c("cov","corcov"),
                          nameToSave = NULL,
                          completePrint=FALSE,
                          buildMethod=c("reg","lassoSS","lassoSSCl"),
                          print=TRUE,
                          cleanCall = FALSE,
                          buildOptions = c("standard"),
                          Rsmlx=c("","2","Step"),
                          numFiles=100){
  if(cleanCall){
    if(print){cat("===========================================================================\n")}
    source("Fun/cleanFile.R")
    for(Rv in Rsmlx){
      cleanFile(covariateSize,covariateType,completePrint,buildMethod,print,buildOptions,Rv)
    }
  }

  for(Rv in Rsmlx){
    for(opt in buildOptions){
      printProj = TRUE
      for(cov in covariateSize){
        for(type in covariateType){
          for(meth in buildMethod){
            printBuildInfo = TRUE
            if(dir.exists(paste0("Results/Results",Rv)) && dir.exists(paste0("Results/Results",Rv,"/Results_",meth))){

              # Folder and paths
              Folder =paste0("Results/Results",Rv,"/")
              resultsFolder = paste0(Folder,"Results_",meth,"/",cov,type)
              pathToSave=paste0(Folder,"Results_",meth,"/finalResults/")
              if(!dir.exists(pathToSave)){dir.create(pathToSave)}

              if(is.null(nameToSave)){
                nameToGiveComp=paste0(cov,type,"comp.RData")
                nameToGive=paste0(cov,type,".RData")
              }else{
                nameToGiveComp=paste0(nameToGive,"comp.RData")
                nameToGive =paste0(nameToGive,".RData")
              }

              # Get Results
              Models = list()
              computationTime = numeric(0)
              iterationCount = numeric(0)
              missingFile = numeric(0)

              for(i in 1:numFiles){
                if(file.exists(paste0(resultsFolder,"/buildResults_",i,".RData"))){
                  load(paste0(resultsFolder,"/buildResults_",i,".RData"))
                  load(paste0(resultsFolder,"/compResults_",i,".RData"))

                  name = paste0("estim_",i)
                  Models <- append(Models,list(model))
                  computationTime <- c(computationTime,time[[1]])
                  iterationCount <- c(iterationCount,iter)

                  names(Models)[length(Models)] <-  names(computationTime)[length(computationTime)] <- names(iterationCount)[length(iterationCount)] <- name
                }else{missingFile = c(missingFile,i)}
              }

              # Print information and save data
              if(length(missingFile)!=numFiles){
                save(Models,file=paste0(pathToSave,nameToGive))
                save(computationTime,iterationCount,file=paste0(pathToSave,"/",nameToGiveComp))

                if(print){
                  if( (length(missingFile)==0 && completePrint) || length(missingFile)!=0){
                    if(printProj){
                      cat("===========================================================================\n")
                      if(opt=="standard"){
                        cat("       - - - - - - - - -  Standard Model Building  - - - - - - - - -       \n")
                      }else if(opt=="Mstar"){
                        cat("       - - - - - - - - - Model Building from M*  - - - - - - - - -       \n")
                      }else{
                        cat(paste0("       - - - - - - - - - Model Building : ",opt,"  - - - - - - - - -       \n"))
                      }
                      printProj=FALSE
                    }
                    if(printBuildInfo){
                      cat(paste0("▶ For project with ",cov))
                      if(type=="corcov"){cat(" correlated covariates, ")}else if(type=="cov"){cat(" covariates,")}
                      cat(" built with ")
                      if(meth == "reg"){
                        cat("stepAIC regression")
                      }else if(meth=="lasso"){
                        cat("lasso regression")
                      }else if(meth=="lassoSS"){
                        cat("lasso with stability selection")
                      }else if(meth=="lassoSSCl"){
                        cat("lasso with stability selection\nand preliminary clusterisation step")
                      }
                      cat(paste0(" ",Rv))
                      cat(" :\n")
                      printBuildInfo = FALSE
                    }
                    if(length(missingFile)!=0){
                      cat(paste0("      • built model missing are : ",paste0(missingFile,collapse=","),".\n"))
                    }else{
                      cat(paste0("      • no built model missing !\n"))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
