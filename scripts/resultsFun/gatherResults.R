gatherResults <- function(project=c("Pasin","GaussianPasin","Naveau"),
                          buildMethod="all",
                          numFiles=100,
                          completePrint=FALSE,
                          print=TRUE){
  
  bM=identical(buildMethod,"all")
  miss = setNames(rep(list(c()),length(project)),project)
  
  for(proj in project){
    if(bM){
      buildMethod <- setdiff(stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",proj),recursive = F),paste0("outputs/buildingResults/simulation/Results",proj,"/Results_")),paste0("outputs/buildingResults/simulation/Results",proj,"/gatheredResults"))
    }
    cat("For ",proj," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
    printProj = TRUE
    
    for(meth in buildMethod){
      if(meth=="SAEMVS"){
          pathToSave ="outputs/buildingResults/simulation/ResultsNaveau/gatheredResults/SAEMVS_finalResults.RData"
          load(paste0("outputs/buildingResults/simulation/ResultsNaveau/gatheredResults/",setdiff(buildMethod,"SAEMVS")[1],"_finalResults.RData"))
          model = Models[[1]]
          computationTime = rep(30*60,numFiles)
          iterationCount = rep(500,numFiles)
          missingFile = c()
          Models <- list() 
          for(i in 1:numFiles){
            if(file.exists(paste0("outputs/buildingResults/simulation/ResultsNaveau/Results_SAEMVS/buildResults_",i,".RData"))){
              load(paste0("outputs/buildingResults/simulation/ResultsNaveau/Results_SAEMVS/buildResults_",i,".RData"))
              
              covariateModel <- list(D=setNames(rep(FALSE,500),paste0("C",1:500)),
                                     V=setNames(rep(FALSE,500),paste0("C",1:500)),
                                     phi1=setNames(rep(FALSE,500),paste0("C",1:500)),
                                     phi2=setNames(rep(FALSE,500),paste0("C",1:500)))
              covariateModel$phi1[res$model_select1] <- TRUE
              covariateModel$phi2[res$model_select2] <- TRUE
              
              formula = paste0(  "log(D) = log(D_pop)\nlog(V) = log(V_pop)\nphi1 = phi1_pop + ",paste0(paste0("beta_phi1_C",res$model_select1,"*C",res$model_select1),collapse=" + ")," + eta_phi1\nphi2 = phi2_pop + ",paste0(paste0("beta_phi2_C",res$model_select2,"*C",res$model_select2),collapse=" + ")," + eta_phi2\nCorrelations\n\tID : {phi1, phi2}\n")
              
              
              
              model$formula <- formula
              model$covariateModel <- covariateModel
              
              Models <- append(Models,list(model))
              
              rm(res)
            }else{missingFile = c(missingFile,i)
            miss[[proj]] = c(miss[[proj]],i)}
          }
          names(Models) <- paste0("estim_",1:numFiles)
          save(Models, computationTime,iterationCount,file=pathToSave)
      }else{
        nameToGive=paste0(meth,"_finalResults.RData")
        
        printBuildInfo = TRUE
        if(dir.exists(paste0("outputs/buildingResults/simulation/Results",proj)) && dir.exists(paste0("outputs/buildingResults/simulation/Results",proj,"/Results_",meth))){
          
          # Folder and paths
          Folder = paste0("outputs/buildingResults/simulation/Results",proj,"/")
          resultsFolder = paste0(Folder,"Results_",meth)
          pathToSave = paste0("outputs/buildingResults/simulation/Results",proj,"/gatheredResults/")
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
            miss[[proj]] = c(miss[[proj]],i)}
          }
          
          # Print information and save data
          if(length(missingFile)!=numFiles){
            save(Models,computationTime,iterationCount,file=paste0(pathToSave,nameToGive))
            
            if(print){
              if( (length(missingFile)==0 && completePrint) || length(missingFile)!=0){
                if(printProj){
                  cat("===========================================================================\n")
                  cat(paste0("       - - - - - - - - -  ",proj," Model Simuations  - - - - - - - - -       \n"))
                  printProj=FALSE
                }
                if(printBuildInfo){
                  cat(paste0("▶ For 200 correlated covariates, "))
                  cat(" built with ", meth)
                  cat(" :\n")
                  printBuildInfo = FALSE
                }
                if(length(missingFile)!=0){
                  cat(paste0("      • built model missing are : ",paste0(sort(unique(missingFile)),collapse=","),".\n"))
                }else{
                  cat(paste0("      • no built model missing !\n"))
                }
              }
            }
          }
        }
      }
      
    }
    if(length(miss[[proj]])!=0){
      cat(paste0("\n\nSUMMARY : missing files ",paste0(sort(unique(miss[[proj]])),collapse=","),".\n\n"))
    }
  }
}
