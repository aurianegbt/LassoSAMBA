suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))

getResults <- function(project=c("Pasin","PK"),
                       buildMethod="all",
                       covariateSize=200,
                       files="all"){
  
  bM = identical(buildMethod,"all")
  
  for(proj in project){
    if(bM){
      buildMethod <- stringr::str_remove_all(list.dirs(paste0("Results/Results",proj),recursive = F),paste0("Results/Results",proj,"/Results_"))
      
      buildMethod = c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassonoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit")[which( c("reg","lasso","elasticnet","lassoSS","elasticnetSS","lassoSSCrit","elasticnetSSCrit",buildMethod[stringr::str_detect(buildMethod,"regPEN")],"regnoCov0","lassonoCov0","elasticnetnoCov0","lassoSSREP","elasticnetSSREP","lassoSSCritREP","elasticnetSSCritREP","rlasso","relasticnet","rsharp","sharp","rlassoCrit","relasticnetCrit") %in% buildMethod)]
    }
    cat("For ",proj," project, build method are ",paste0(buildMethod,collapse = ", "),"\n")
    
    
    source(paste0("Files/Files",proj,"/H1.all.R"))
    
    
    if(file.exists(paste0("Save/orderList",proj,".RData"))){
      load(paste0("Save/orderList",proj,".RData"))
    }else{
      orderList = list()
      
      for(type in c("cov","corcov")){
        sim = max(list.dirs(paste0("Files/Files",proj,"/"),recursive=FALSE) %>%
                    str_remove_all(paste0("Files/Files",proj,"/")) %>%
                    str_remove_all("corcov") %>%
                    str_remove_all("cov") %>%
                    setdiff("keep"))
        covTable = read.csv(paste0("Files/Files",proj,"/",sim,type,"/covTable/covTable_1.txt"))
        covnames = colnames(covTable)[-1]
        
        orderList = append(orderList,list(covnames))
        names(orderList)[length(orderList)] <- type
      }
      
      save(orderList,file=paste0("Save/orderList",proj,".RData"))
    }
    if(identical(files,"all")){files=1:100}
    if(identical(files,"all.computed")){
      files.list=list()
      for(meth in buildMethod){
        for(sim in covariateSize){
          for(type in c("cov","corcov")){
            pathToRes =paste0("Results/Results",proj,"/Results_",meth,"/finalResults/",sim,type,".RData")
            
            if(file.exists(pathToRes)){
              load(pathToRes)
              files.list = append(files.list,list(as.numeric(stringr::str_remove(names(Models),"estim_"))))
              names(files.list)[length(files.list)] <- paste0(meth,sim,cov)
            }
          }
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
      for(sim in covariateSize){
        for(type in c("cov","corcov")){
          pathToRes =paste0("Results/Results",proj,"/Results_",meth,"/finalResults/")
          pathToCompRes = paste0(pathToRes,sim,type,"comp.RData")
          pathToRes = paste0(pathToRes,sim,type,".RData")
          if(file.exists(pathToRes)){
            load(pathToRes) #Models
            load(pathToCompRes) #computationTime, iterationCount
            
            n = length(Models)
            to.consider = intersect(paste0("estim_",files),names(Models))
            ind.consider = which(names(Models) %in% to.consider)
            
            
            TotalNumberofModel=c(TotalNumberofModel,length(to.consider))
            names(TotalNumberofModel)[length(TotalNumberofModel)] <- paste0(meth,sim,type)
            
            selectionPar = list()
            
            for (j in ind.consider){
              Model = Models[[j]]
              num = stringr::str_replace(names(Models)[j],"estim_","")
              
              for(p in param){
                if(length(names(which(Model$covariateModel[[p]])))!=0){
                  CovariateModelSelection <- rbind(CovariateModelSelection,
                                                   data.frame(Covariate = names(which(Model$covariateModel[[p]])),
                                                              Parameter = t.param[[p]],
                                                              Model =num,ProjectNumber=paste0(sim," covariates"),
                                                              TypeOfSim = type,
                                                              Method = meth))
                }
              }
              
              selectionPar = lapply(Model$covariateModel,function(x){names(which(x))})
              selection = Reduce(union,selectionPar)
              
              H1 = Reduce(union,H1.all)
              H0 = setdiff(orderList[[type]][1:as.numeric(sim)],H1)
              
              
              errorStats <- rbind(errorStats,data.frame(Model = num,
                                                        ProjectNumber=paste0(sim," covariates"),
                                                        TypeOfSim = type,
                                                        Method = meth,
                                                        FP = FP(selection,H0,H1),
                                                        TP = TP(selection,H0,H1),
                                                        FDR = FDR(selection,H0,H1),
                                                        FN = FN(selection,H0,H1),
                                                        TN = TN(selection,H0,H1),
                                                        FNR = FNR(selection,H0,H1)))
              
              H0.all = lapply(H1.all,FUN = function(x){setdiff(orderList[[type]][1:as.numeric(sim)],x)})
              
              aux = lapply(names(H1.all),FUN = function(x){data.frame(Model = num,
                                                                      ProjectNumber=paste0(sim," covariates"),
                                                                      TypeOfSim = type,
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
                                                                 ProjectNumber=paste0(sim," covariates"),
                                                                 TypeOfSim = type,
                                                                 Method = meth))
          }
        }
      }
    }    
        
    resultCovariate <- data.frame()
    resultCovariatePar <- data.frame()
    resultModelPar <- data.frame()
    resultModel <- data.frame()
    
    for(type in c("cov","corcov")){
      covnames = orderList[[type]]
      nbCov = length(covnames)
      
      for(i in 1:nbCov){
        covariate = covnames[i]
        df = CovariateModelSelection[CovariateModelSelection$Covariate==covariate,]
        for(meth in buildMethod){
          for(sim in covariateSize){
            if(paste0(meth,sim,type) %in% names(TotalNumberofModel)){
              NumberofModel = TotalNumberofModel[[paste0(meth,sim,type)]]
              aux = df[df$ProjectNumber==paste0(sim," covariates") &
                         df$TypeOfSim==type &
                         df$Method==meth,]
              SelectedInDistinctModel = length(unique(aux$Model))
              resultCovariate <- rbind(resultCovariate,data.frame(Covariate =covariate,
                                                                  ProjectNumber = paste0(sim," covariates"),
                                                                  TypeOfSim = type,
                                                                  Method = meth,
                                                                  NumberofModel=NumberofModel,
                                                                  SelectedInDistinctModel=SelectedInDistinctModel/NumberofModel))
              ## Par
              for(par in unname(t.param)){
                aux = df[df$ProjectNumber==paste0(sim," covariates") &
                           df$TypeOfSim==type &
                           df$Method==meth & df$Parameter==par,]
                ProportionSelected = nrow(aux) / NumberofModel
                
                resultCovariatePar <- rbind(resultCovariatePar,data.frame(Covariate =covariate,
                                                                          Parameter = par,
                                                                          ProjectNumber = paste0(sim," covariates"),
                                                                          TypeOfSim = type,
                                                                          Method = meth,
                                                                          NumberofModel=NumberofModel,
                                                                          ProportionSelected=ProportionSelected))
              }
            }
          }
        }
      }
      
      for(meth in buildMethod){
        for(sim in covariateSize){
          if(paste0(meth,sim,type) %in% names(TotalNumberofModel)){
            
            NoFNModel = length(unique(errorStats[errorStats$FN==0 &
                                                   errorStats$ProjectNumber==paste0(sim," covariates") &
                                                   errorStats$TypeOfSim==type &
                                                   errorStats$Method==meth,"Model"]))
            TrueModel = length(unique(errorStats[errorStats$FN==0 &
                                                   errorStats$FP==0 &
                                                   errorStats$ProjectNumber==paste0(sim," covariates") &
                                                   errorStats$TypeOfSim==type &
                                                   errorStats$Method==meth,"Model"]))
            
            resultModel <- rbind(resultModel,data.frame(ProjectNumber = paste0(sim," covariates"),
                                                        TypeOfSim = type,
                                                        Method = meth,
                                                        TrueModel = TrueModel/TotalNumberofModel[[paste0(meth,sim,type)]],
                                                        NoFNModel = NoFNModel/TotalNumberofModel[[paste0(meth,sim,type)]]))
            
            errorStatsParTrans = suppressMessages(errorStatsPar %>%
              group_by(Model,ProjectNumber,TypeOfSim,Method) %>%
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
                                                   errorStatsParTrans$ProjectNumber==paste0(sim," covariates") &
                                                   errorStatsParTrans$TypeOfSim==type &
                                                   errorStatsParTrans$Method==meth,"Model"]))
            TrueModel = length(unique(errorStatsParTrans[errorStatsParTrans$FN==0 &
                                                   errorStatsParTrans$FP==0 &
                                                   errorStatsParTrans$ProjectNumber==paste0(sim," covariates") &
                                                   errorStatsParTrans$TypeOfSim==type &
                                                   errorStatsParTrans$Method==meth,"Model"]))
            
            resultModelPar <- rbind(resultModelPar,data.frame(ProjectNumber = paste0(sim," covariates"),
                                                        TypeOfSim = type,
                                                        Method = meth,
                                                        TrueModel = TrueModel/TotalNumberofModel[[paste0(meth,sim,type)]],
                                                        NoFNModel = NoFNModel/TotalNumberofModel[[paste0(meth,sim,type)]]))
            
            
          }
        }
      }
    }
    
    
    save(computationStats,errorStats,errorStatsPar,resultCovariate,orderList,CovariateModelSelection,resultCovariatePar,resultModel,resultModelPar,file=paste0("Save/BuildResults_",proj,".RData"))
  }
} 
