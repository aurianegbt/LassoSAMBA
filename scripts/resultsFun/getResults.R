suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))

getResults <- function(project=c("Pasin","GaussianPasin"),
                       buildMethod="all",
                       exclude=c("sharp","sharpnoCov0"),
                       files="all"){
  
  bM = identical(buildMethod,"all")
  
  for(proj in project){
    source(paste0("data/simulationFiles/Files",proj,"/H1.all.R"))
    if(bM){
      buildMethod <- stringr::str_remove_all(list.dirs(paste0("outputs/buildingResults/simulation/Results",proj),recursive = F),paste0("outputs/buildingResults/simulation/Results",proj,"/Results_"))
      buildMethod = c("reg",
                      buildMethod[stringr::str_detect(buildMethod,"lassonoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"lassonoCov0FDP"))])[which(c("reg", buildMethod[stringr::str_detect(buildMethod,"lassonoCov0FDP") & grepl("^[0-9]+$", stringr::str_remove_all(buildMethod,"lassonoCov0FDP"))]) %in% buildMethod)]
      if(!is.null(exclude)){
        buildMethod = setdiff(buildMethod,exclude)
      }
      
      
    }
    cat("For ",proj," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
    
    
    if(file.exists(paste0("data/simulationFiles/orderList",proj,".RData"))){
      load(paste0("data/simulationFiles/orderList",proj,".RData"))
    }else{
      load(paste0("data/simulationFiles/Files",proj,"/headerTypes.RData"))
      
      orderList = colnames(read.csv(paste0("data/simulationFiles/Files",proj,"/simulation/simulation_1.txt")))[which(stringr::str_detect(headerTypes,"cov"))]
      
      save(orderList,file=paste0("data/simulationFiles/orderList",proj,".RData"))
    }
    if(identical(files,"all")){files=1:100}
    if(identical(files,"all.computed")){
      files.list=list()
      for(meth in buildMethod){
        pathToRes = paste0("outputs/buildingResults/simulation/Results",proj,"/finalResults.RData")
        
        if(file.exists(pathToRes)){
          load(pathToRes)
          files.list = append(files.list,list(as.numeric(stringr::str_remove(names(Models),"estim_"))))
          names(files.list)[length(files.list)] <- paste0(meth)
        }
      }
      files = 1:100
      for(l in files.list){
        files = intersect(files,l)
      }
      if(length(files)<100){
        warning("Less than 100 files !")
      }
    }
    
    TP <- function(sel,H0,H1){length(unique(intersect(sel,H1)))}
    TN <- function(sel,H0,H1){length(unique(setdiff(H0,sel)))}
    FP <- function(sel,H0,H1){length(unique(intersect(sel,H0)))}
    FN <- function(sel,H0,H1){length(unique(setdiff(H1,sel)))}
    
    FDR <- function(sel,H0,H1){round(FP(sel,H0,H1)/(max(FP(sel,H0,H1)+TP(sel,H0,H1),1)),digits=3)}
    FNR <- function(sel,H0,H1){round(FN(sel,H0,H1)/(max(FN(sel,H0,H1)+TN(sel,H0,H1),1)),digits=3)}
    
        
    CovariateModelSelection <- data.frame()
    TotalNumberofModel=c()
    computationStats <- data.frame()
    errorStats <- data.frame()
    errorStatsPar <- data.frame()
    
    for(meth in buildMethod){
      pathToRes = paste0("outputs/buildingResults/simulation/Results",proj,"/gatheredResults/",meth,"_finalResults.RData")
      if(file.exists(pathToRes)){
        load(pathToRes) #Models, computationTime, iterationCount
        
        n = length(Models)
        to.consider = intersect(paste0("estim_",files),names(Models))
        ind.consider = which(names(Models) %in% to.consider)
        
        
        TotalNumberofModel=c(TotalNumberofModel,length(to.consider))
        names(TotalNumberofModel)[length(TotalNumberofModel)] <- paste0(meth)
        
        selectionPar = list()
        
        for (j in ind.consider){
          Model = Models[[j]]
          num = stringr::str_replace(names(Models)[j],"estim_","")
          
          for(p in param){
            if(length(names(which(Model$covariateModel[[p]])))!=0){
              CovariateModelSelection <- rbind(CovariateModelSelection,
                                               data.frame(Model = num,
                                                          Covariate = names(which(Model$covariateModel[[p]])),
                                                          Parameter = t.param[[p]],
                                                          Method = meth))
            }
          }
          
          selectionPar = lapply(Model$covariateModel,function(x){names(which(x))})
          selection = Reduce(union,selectionPar)
          
          H1 = Reduce(union,H1.all)
          H0 = setdiff(orderList,H1)
          
          
          errorStats <- rbind(errorStats,data.frame(Model = num,
                                                    Method = meth,
                                                    FP = FP(selection,H0,H1),
                                                    TP = TP(selection,H0,H1),
                                                    FDR = FDR(selection,H0,H1),
                                                    FN = FN(selection,H0,H1),
                                                    TN = TN(selection,H0,H1),
                                                    FNR = FNR(selection,H0,H1)))
          
          H0.all = lapply(H1.all,FUN = function(x){setdiff(orderList,x)})
          
          aux = lapply(names(H1.all),FUN = function(x){data.frame(Model = num,
                                                                  Method = meth,
                                                                  Parameter = x,
                                                                  FP = FP(selectionPar[[x]],H0.all[[x]],H1.all[[x]]),
                                                                  TP = TP(selectionPar[[x]],H0.all[[x]],H1.all[[x]]),
                                                                  FDR = FDR(selectionPar[[x]],H0.all[[x]],H1.all[[x]]),
                                                                  FN = FN(selectionPar[[x]],H0.all[[x]],H1.all[[x]]),
                                                                  TN = TN(selectionPar[[x]],H0.all[[x]],H1.all[[x]]),
                                                                  FNR = FNR(selectionPar[[x]],H0.all[[x]],H1.all[[x]]))})
          
          errorStatsPar <- do.call("rbind",append(list(errorStatsPar),aux))
        }
        computationStats <-rbind(computationStats,data.frame(time = computationTime[ind.consider],
                                                             iteration = iterationCount[ind.consider],
                                                             Model = stringr::str_replace(names(computationTime),"estim_","")[ind.consider],
                                                             Method = meth))
      }
    }
        
    resultCovariate <- data.frame()
    resultCovariatePar <- data.frame()
    resultModelPar <- data.frame()
    resultModel <- data.frame()
    
    for(i in 1:length(orderList)){
      covariate = orderList[i]
      df = CovariateModelSelection[CovariateModelSelection$Covariate==covariate,]
      for(meth in buildMethod){
        if(meth %in% names(TotalNumberofModel)){
          NumberofModel = TotalNumberofModel[[meth]]
          aux = df[df$Method==meth,]
          SelectedInDistinctModel = length(unique(aux$Model))
          resultCovariate <- rbind(resultCovariate,data.frame(Covariate = covariate,
                                                              Method = meth,
                                                              NumberofModel=NumberofModel,
                                                              SelectedInDistinctModel=SelectedInDistinctModel/NumberofModel))
          ## Par
          for(par in unname(t.param)){
            aux = df[df$Method==meth & df$Parameter==par,]
            ProportionSelected = nrow(aux) / NumberofModel
            
            resultCovariatePar <- rbind(resultCovariatePar,data.frame(Covariate = covariate,
                                                                      Parameter = par,
                                                                      Method = meth,
                                                                      NumberofModel=NumberofModel,
                                                                      ProportionSelected=ProportionSelected))
          }
        }
      }
    }
      
    for(meth in buildMethod){
      if(meth %in% names(TotalNumberofModel)){
        
        NoFNModel = length(unique(errorStats[errorStats$FN==0 &
                                               errorStats$Method==meth,"Model"]))
        TrueModel = length(unique(errorStats[errorStats$FN==0 &
                                               errorStats$FP==0 &
                                               errorStats$Method==meth,"Model"]))
        
        resultModel <- rbind(resultModel,data.frame(Method = meth,
                                                    TrueModel = TrueModel/TotalNumberofModel[[meth]],
                                                    NoFNModel = NoFNModel/TotalNumberofModel[[meth]]))
        
        errorStatsParTrans = suppressMessages(errorStatsPar %>%
                                                group_by(Model,Method) %>%
                                                summarise(
                                                  FP = sum(FP),
                                                  TP = sum(TP), 
                                                  FN = sum(FN),
                                                  TN = sum(TN)
                                                ) %>% 
                                                mutate(FDR = (FP/max(FP+TP,1)),.after = "TP") %>%
                                                mutate(FNR = (FN/max(FN+TN,1)),.after="TN") %>%
                                                as.data.frame())
        
        NoFNModel = length(unique(errorStatsParTrans[errorStatsParTrans$FN==0 &
                                                       errorStatsParTrans$Method==meth,"Model"]))
        TrueModel = length(unique(errorStatsParTrans[errorStatsParTrans$FN==0 &
                                                       errorStatsParTrans$FP==0 &
                                                       errorStatsParTrans$Method==meth,"Model"]))
        
        resultModelPar <- rbind(resultModelPar,data.frame(Method = meth,
                                                          TrueModel = TrueModel/TotalNumberofModel[[meth]],
                                                          NoFNModel = NoFNModel/TotalNumberofModel[[meth]]))
        
        
      }
    }
    
    save(computationStats,errorStats,errorStatsPar,resultCovariate,orderList,CovariateModelSelection,resultCovariatePar,resultModel,resultModelPar,file=paste0("outputs/finalResults/BuildResults_",proj,".RData"))
  }
} 
