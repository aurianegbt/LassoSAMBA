gatherResults <- function(project=c("Pasin","GaussianPasin"),
                          buildMethod="all",
                          exclude=c("sharp","sharpnoCov0"),
                          numFiles=100,
                          completePrint=FALSE,
                          print=TRUE,
                          cleanCall = TRUE){
  
  bM=identical(buildMethod,"all")
  miss = setNames(rep(list(c()),length(project)),project)
  
  for(proj in project){
    if(bM){
      buildMethod <- stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",proj),recursive = F),paste0("outputs/buildingResults/simulation/Results",proj,"/Results_"))
      buildMethod = c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharpnoCov0"))],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharpnoCov0FDP"))],buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharp"))],"regnoCov0","lassoSSnoCov0","lassoSSCritnoCov0","sharpnoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit")[which( c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharpnoCov0"))],buildMethod[stringr::str_detect(buildMethod,"sharpnoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharpnoCov0FDP"))],buildMethod[stringr::str_detect(buildMethod,"sharp") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"sharp"))],"regnoCov0","lassoSSnoCov0","lassoSSCritnoCov0","sharpnoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit") %in% buildMethod)]
      if(!is.null(exclude)){
        buildMethod = setdiff(buildMethod,exclude)
      }
    }
    cat("For ",proj," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
    
    
    
    if(cleanCall){
      source("scripts/resultsFun/cleanFile.R")
      cleanFile(proj,buildMethod,completePrint,print)
    }
    printProj = TRUE
    
    for(meth in buildMethod){
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
    if(length(miss[[proj]])!=0){
      cat(paste0("\n\nSUMMARY : missing files ",paste0(sort(unique(miss[[proj]])),collapse=","),".\n\n"))
    }
  }
}
