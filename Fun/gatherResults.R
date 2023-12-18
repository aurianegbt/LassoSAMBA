gatherResults <- function(project=c("Pasin","Warfarine"),
                          covariateSize=c(10,50,200,500),
                          buildMethod=c("reg","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit"),
                          numFiles=100,
                          completePrint=FALSE,
                          print=TRUE,
                          cleanCall = TRUE){
  if(cleanCall){
    if(print){cat("===========================================================================\n")}
    source("Fun/cleanFile.R")
      cleanFile(project,covariateSize,buildMethod,completePrint,print)
  }
  
  for(proj in project){
    miss = setNames(rep(list(c()),length(project)),project)
    printProj = TRUE
      for(cov in covariateSize){
        for(type in c("cov","corcov")){
          
          nameToGiveComp=paste0(cov,type,"comp.RData")
          nameToGive=paste0(cov,type,".RData")
          
          for(meth in buildMethod){
            printBuildInfo = TRUE
            if(dir.exists(paste0("Results/Results",proj)) && dir.exists(paste0("Results/Results",proj,"/Results_",meth))){

              # Folder and paths
              Folder =paste0("Results/Results",proj,"/")
              resultsFolder = paste0(Folder,"Results_",meth,"/",cov,type)
              pathToSave=paste0(Folder,"Results_",meth,"/finalResults/")
              if(!dir.exists(pathToSave)){dir.create(pathToSave)}


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
                }else{missingFile = c(missingFile,i)
                miss[[proj]] <- sort(unique(c(miss[[proj]],i)))
                }
              }

              # Print information and save data
              if(length(missingFile)!=numFiles){
                save(Models,file=paste0(pathToSave,nameToGive))
                save(computationTime,iterationCount,file=paste0(pathToSave,"/",nameToGiveComp))

                if(print){
                  if( (length(missingFile)==0 && completePrint) || length(missingFile)!=0){
                    if(printProj){
                      cat("===========================================================================\n")
                      cat(paste0("       - - - - - - - - -  ",proj," Model Simuations  - - - - - - - - -       \n"))
                      printProj=FALSE
                    }
                    if(printBuildInfo){
                      cat(paste0("▶ For ",cov))
                      if(type=="corcov"){cat(" correlated covariates, ")}else if(type=="cov"){cat(" covariates,")}
                      cat(" built with ")
                      if(meth == "reg"){
                        cat("stepAIC regression")
                      }else if(stringr::str_detect(meth,"lasso")){
                        cat("lasso")
                      }else if(stringr::str_detect(meth,"elasticnet")){
                        cat("elastic net")
                      }
                      if(stringr::str_detect(meth,"SS")){
                        cat(" with stability selection")
                      }
                      if(stringr::str_detect(meth,"Crit")){
                        cat("\nwith multiple thresholds")
                      }
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
    if(length(miss[[proj]])!=0){
      cat(paste0("\n\nSUMMARY : missing files ",paste0(miss[[proj]],collapse=","),".\n\n"))
    }
    }
  }
