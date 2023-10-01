getResults <- function(buildMethod=c("reg","lassoSS","lassoSSCl"),
                       covariateSize=c(100,200),
                       covariateType=c("cov","corcov"),
                       files="all.computed",
                       buildOptions="standard",
                       Rsmlx=c("","2","Step")){
  getResults.data(buildMethod,covariateSize,covariateType,files,buildOptions,Rsmlx)
  getResults.pardata(buildMethod,covariateSize,covariateType,files,buildOptions,Rsmlx)
}


getResults.data <- function(buildMethod=c("reg","lassoSS","lassoSSCl"),
                              covariateSize=c(100,200),
                              covariateType=c("cov","corcov"),
                              files="all.computed",
                              buildOptions="standard",
                              Rsmlx=c("","2","Step")){
  if(file.exists("Save/orderList.RData")){
    load("Save/orderList.RData")
  }else{
    orderList = list()

    for(type in covariateType){
      sim = max(covariateSize)
      covTable = read.csv(paste0("Files/",sim,stringr::str_remove(type,"2"),"/covTable/covTable_1.txt"))
      covnames = colnames(covTable)[-1]

      load(paste0("Files/",sim,stringr::str_remove(type,"2"),"/headerTypes.RData"))
      contCov = which(headerTypes[-c(1:3)]=="contcov")

      covnames[contCov] <- paste0("c",covnames[contCov])

      orderList = append(orderList,list(covnames))
      names(orderList)[length(orderList)] <- type
    }

    save(orderList,file="Save/orderList.RData")
  }
  if(identical(files,"all")){files=1:100}
  if(identical(files,"all.computed")){
    files.list=list()
    for(Rv in Rsmlx){
      for(opt in buildOptions){
        for(meth in buildMethod){
          for(sim in covariateSize){
            for(type in covariateType){
              pathToRes =paste0("Results/Results",Rv,"/","Results_",meth,"/finalResults/",sim,type,".RData")

              if(file.exists(pathToRes)){
                load(pathToRes)
                files.list = append(files.list,list(as.numeric(stringr::str_remove(names(Models),"estim_"))))
                names(files.list)[length(files.list)] <- paste0(meth,opt,type)
              }
            }
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

  for(Rv in Rsmlx){
    for(opt in buildOptions){

      CovariateModelSelection <- data.frame()
      TotalNumberofModel=c()
      computationStats <- data.frame()
      errorStats <- data.frame()
      errorStatsPar <- data.frame()


      for(meth in buildMethod){
        for(sim in covariateSize){
          for(type in covariateType){
            pathToRes =paste0("Results/Results",Rv,"/")
            pathToRes = paste0(pathToRes,"Results_",meth,"/finalResults/")
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


              for (j in ind.consider){
                Model = Models[[j]]
                num = stringr::str_replace(names(Models)[j],"estim_","")


                if(length(names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]))!=0){
                  auxS = data.frame(Covariate = names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]),
                                    Parameter = latex2exp::TeX(r"($\varphi_S$)",output="character"),
                                    Model =num,ProjectNumber=paste0(sim," covariates"),
                                    TypeOfSim = type,
                                    Method = meth)
                }else{auxS = data.frame()}
                if(length(names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]))!=0){
                  auxL = data.frame(Covariate = names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]),
                                    Parameter = latex2exp::TeX(r"($\varphi_L$)",output="character"),
                                    Model=num,
                                    ProjectNumber=paste0(sim," covariates"),
                                    TypeOfSim = type,
                                    Method = meth)
                }else{ auxL = data.frame()}
                if(length(names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]))!=0){
                  auxD =data.frame(Covariate = names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]),
                                   Parameter = latex2exp::TeX(r"($\delta_{Ab}$)",output="character"),
                                   Model=num,
                                   ProjectNumber=paste0(sim," covariates"),
                                   TypeOfSim = type,
                                   Method = meth)
                }else{auxD = data.frame()}
                CovariateModelSelection <- rbind(CovariateModelSelection,auxS,auxL,auxD)


                ## phi_S
                selectionS = names(Model$covariateModel$phi_S[Model$covariateModel$phi_S])
                selectionL = names(Model$covariateModel$phi_L[Model$covariateModel$phi_L])
                selectionD = names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB])

                selection = union(selectionS,union(selectionL,selectionD))
                H1 = c("cAGE","SEX","RACE")
                H0 = setdiff(orderList[[stringr::str_remove(type,"2")]][1:as.numeric(sim)],H1)


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

                H1 = list(phi_S="cAGE",phi_L="RACE",delta_AB="SEX")
                H0 = lapply(H1,FUN = function(x){setdiff(orderList[[stringr::str_remove(type,"2")]][1:as.numeric(sim)],x)})

                aux = lapply(list("phi_S","phi_L","delta_AB"),FUN = function(x){data.frame(Model = num,
                                                              ProjectNumber=paste0(sim," covariates"),
                                                              TypeOfSim = type,
                                                              Method = meth,
                                                              Parameter = x,
                                                              FP = FP(selection[[x]],H0[[x]],H1[[x]]),
                                                              TP = TP(selection[[x]],H0[[x]],H1[[x]]),
                                                              FDR = FDR(selection[[x]],H0[[x]],H1[[x]]),
                                                              FN = FN(selection[[x]],H0[[x]],H1[[x]]),
                                                              TN = TN(selection[[x]],H0[[x]],H1[[x]]),
                                                             FNR = FNR(selection[[x]],H0[[x]],H1[[x]]))})

                errorStatsPar <- rbind(errorStatsPar,aux[[1]],aux[[2]],aux[[3]])


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
      resultModel <- data.frame()
      covLink = c(cAGE = latex2exp::TeX(r"($\varphi_S$)",output="character"),
                  RACE=latex2exp::TeX(r"($\varphi_L$)",output="character"),
                  SEX=latex2exp::TeX(r"($\delta_{Ab}$)",output="character"))

      for(type in covariateType){
        covnames = orderList[[stringr::str_remove(type,"2")]]
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
                for(par in c("varphi[S]","varphi[L]","delta[Ab]")){
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
              aux = CovariateModelSelection[CovariateModelSelection$Covariate %in% c("cAGE","RACE","SEX") & CovariateModelSelection$ProjectNumber==paste0(sim," covariates") & CovariateModelSelection$TypeOfSim==type & CovariateModelSelection$Method==meth,]

              NoFNModel = 0
              for(model in unique(aux$Model)){
                if(nrow(aux[aux$Model==model,])==3 && identical(aux[aux$Model==model & aux$Covariate=="cAGE","Parameter"],"varphi[S]") && identical(aux[aux$Model==model & aux$Covariate=="RACE","Parameter"],"varphi[L]") && identical(aux[aux$Model==model & aux$Covariate=="SEX","Parameter"],"delta[Ab]")){
                  NoFNModel = NoFNModel +1
                }
              }
              aux = CovariateModelSelection[CovariateModelSelection$ProjectNumber==paste0(sim," covariates") & CovariateModelSelection$TypeOfSim==type & CovariateModelSelection$Method==meth,]

              TrueModel = 0
              for(model in unique(aux$Model)){
                if(setequal(aux[aux$Model==model,"Covariate"],c("cAGE","SEX","RACE")) && identical(aux[aux$Model==model & aux$Covariate=="cAGE","Parameter"],"varphi[S]") && identical(aux[aux$Model==model & aux$Covariate=="RACE","Parameter"],"varphi[L]") && identical(aux[aux$Model==model & aux$Covariate=="SEX","Parameter"],"delta[Ab]")){
                  TrueModel = TrueModel + 1
                }
              }



              resultModel <- rbind(resultModel,data.frame(ProjectNumber = paste0(sim," covariates"),
                                                          TypeOfSim = type,
                                                          Method = meth,
                                                          TrueModel = TrueModel/TotalNumberofModel[[paste0(meth,sim,type)]],
                                                          NoFNModel = NoFNModel/TotalNumberofModel[[paste0(meth,sim,type)]]))
            }
          }
        }
      }
      resultCovariate[resultCovariate$Method=="reg","Method"] <- "regression"
      resultCovariatePar[resultCovariatePar$Method=="reg","Method"] <- "regression"
      resultModel[resultModel$Method=="reg","Method"] <- "regression"
      CovariateModelSelection[CovariateModelSelection$Method=="reg","Method"] <- "regression"
      computationStats[computationStats$Method=="reg","Method"] <- "regression"
      errorStats[errorStats$Method=="reg","Method"] <- "regression"

      resultCovariate[resultCovariate$Method=="CoVSS","Method"] <- "ClustOfVar"
      resultCovariatePar[resultCovariatePar$Method=="CoVSS","Method"] <- "ClustOfVar"
      resultModel[resultModel$Method=="CoVSS","Method"] <- "ClustOfVar"
      CovariateModelSelection[CovariateModelSelection$Method=="CoVSS","Method"] <- "ClustOfVar"
      computationStats[computationStats$Method=="CoVSS","Method"] <- "ClustOfVar"
      errorStats[errorStats$Method=="CoVSS","Method"] <- "ClustOfVar"


      save(computationStats,errorStats,errorStatsPar,resultCovariate,orderList,CovariateModelSelection,resultCovariatePar,resultModel,file=paste0("Save/BuildResults_",Rv,".RData"))
    }
  }
}

#===============================================================================================================#

getResults.pardata <- function(buildMethod=c("reg","lassoSS","lassoSSCl"),
                               covariateSize=c(100,200),
                               covariateType=c("cov","corcov"),
                               files="all.computed",
                               buildOptions="standard",
                               Rsmlx=c("","2","Step")){
  if(file.exists("Save/orderList.RData")){
    load("Save/orderList.RData")
  }else{
    orderList = list()

    for(type in covariateType){
      sim = max(covariateSize)
      covTable = read.csv(paste0("Files/",sim,stringr::str_remove(type,"2"),"/covTable/covTable_1.txt"))
      covnames = colnames(covTable)[-1]

      load(paste0("Files/",sim,stringr::str_remove(type,"2"),"/headerTypes.RData"))
      contCov = which(headerTypes[-c(1:3)]=="contcov")

      covnames[contCov] <- paste0("c",covnames[contCov])

      orderList = append(orderList,list(covnames))
      names(orderList)[length(orderList)] <- type
    }

    save(orderList,file="Save/orderList.RData")
  }
  if(identical(files,"all")){files=1:100}
  if(identical(files,"all.computed")){
    files.list=list()
    for(Rv in Rsmlx){
      for(opt in buildOptions){
        for(meth in buildMethod){
          for(sim in covariateSize){
            for(type in covariateType){
              pathToRes =paste0("Results/Results",Rv,"/","Results_",meth,"/finalResults/",sim,type,".RData")

              if(file.exists(pathToRes)){
                load(pathToRes)
                files.list = append(files.list,list(as.numeric(stringr::str_remove(names(Models),"estim_"))))
                names(files.list)[length(files.list)] <- paste0(meth,opt,type)
              }
            }
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

  for(Rv in Rsmlx){
    for(opt in buildOptions){
      CovariateModelSelection <- data.frame()
      TotalNumberofModel=c()

      for(meth in buildMethod){
        for(sim in covariateSize){
          for(type in covariateType){
            pathToRes =paste0("Results/Results",Rv,"/")
            pathToRes = paste0(pathToRes,"Results_",meth,"/finalResults/",sim,type,".RData")
            if(file.exists(pathToRes)){
              load(pathToRes) #Models

              n = length(Models)
              to.consider = intersect(paste0("estim_",files),names(Models))
              ind.consider = which(names(Models) %in% to.consider)


              TotalNumberofModel=c(TotalNumberofModel,length(to.consider))
              names(TotalNumberofModel)[length(TotalNumberofModel)] <- paste0(meth,sim,type)


              for (j in ind.consider){
                Model = Models[[j]]
                num = stringr::str_replace(names(Models)[j],"estim_","")

                if(length(names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]))!=0){
                  auxS = data.frame(Covariate = names(Model$covariateModel$phi_S[Model$covariateModel$phi_S]),
                                    Parameter = latex2exp::TeX(r"($\varphi_S$)",output="character"),
                                    Model=num,
                                    ProjectNumber=latex2exp::TeX(paste0(sim," covariates"),output="character"),
                                    TypeOfSim = type,
                                    Method = meth)
                }else{auxS = data.frame()}
                if(length(names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]))!=0){
                  auxL = data.frame(Covariate = names(Model$covariateModel$phi_L[Model$covariateModel$phi_L]),
                                    Parameter = latex2exp::TeX(r"($\varphi_L$)",output="character"),
                                    Model=num,
                                    ProjectNumber=latex2exp::TeX(paste0(sim," covariates"),output="character"),
                                    TypeOfSim = type,
                                    Method = meth)
                }else{ auxL = data.frame()}
                if(length(names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]))!=0){
                  auxD =data.frame(Covariate = names(Model$covariateModel$delta_AB[Model$covariateModel$delta_AB]),
                                   Parameter = latex2exp::TeX(r"($\delta_{Ab}$)",output="character"),
                                   Model=num,
                                   ProjectNumber=latex2exp::TeX(paste0(sim," covariates"),output="character"),
                                   TypeOfSim = type,
                                   Method = meth)
                }else{auxD = data.frame()}
                CovariateModelSelection <- rbind(CovariateModelSelection,auxS,auxL,auxD)
              }
            }
          }
        }
      }


      resultCovariate <- data.frame()
      resultCovariatePar <- data.frame()
      resultModel <- data.frame()
      covLink = c(cAGE = latex2exp::TeX(r"($\varphi_S$)",output="character"),
                  RACE=latex2exp::TeX(r"($\varphi_L$)",output="character"),
                  SEX=latex2exp::TeX(r"($\delta_{Ab}$)",output="character"))

      for(type in covariateType){
        covnames = orderList[[stringr::str_remove(type,"2")]]
        nbCov = length(covnames)

        for(i in 1:nbCov){
          covariate = covnames[i]
          df = CovariateModelSelection[CovariateModelSelection$Covariate==covariate,]
          for(meth in buildMethod){
            for(sim in covariateSize){
              if(paste0(meth,sim,type) %in% names(TotalNumberofModel)){
                NumberofModel = TotalNumberofModel[[paste0(meth,sim,type)]]
                aux = df[df$ProjectNumber==latex2exp::TeX(paste0(sim," covariates"),output="character") &
                           df$TypeOfSim==type &
                           df$Method==meth,]
                TotalSelected = nrow(aux)
                SelectedInDistinctModel = length(unique(aux$Model))
                resultCovariate <- rbind(resultCovariate,
                                         data.frame(Covariate =covariate,
                                                    ProjectNumber = latex2exp::TeX(paste0(sim," covariates"),output="character"),
                                                    TypeOfSim = type,
                                                    Method = meth,
                                                    NumberofModel=NumberofModel,
                                                    SelectedInDistinctModel=SelectedInDistinctModel/NumberofModel))

                ## Par
                NumberofModel = TotalNumberofModel[[paste0(meth,sim,type)]]
                for(par in c("varphi[S]","varphi[L]","delta[Ab]")){
                  aux = df[df$ProjectNumber==latex2exp::TeX(paste0(sim," covariates"),output="character") &
                             df$TypeOfSim==type &
                             df$Method==meth & df$Parameter==par,]
                  ProportionSelected = nrow(aux) / NumberofModel

                  resultCovariatePar <- rbind(resultCovariatePar,data.frame(Covariate =covariate,
                                                                            Parameter = par,
                                                                            ProjectNumber = latex2exp::TeX(paste0(sim," covariates"),output="character"),
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
              aux = CovariateModelSelection[CovariateModelSelection$Covariate %in% c("cAGE","RACE","SEX") & CovariateModelSelection$ProjectNumber==latex2exp::TeX(paste0(sim," covariates"),output="character")  & CovariateModelSelection$TypeOfSim==type & CovariateModelSelection$Method==meth,]

              NoFNModel = 0
              for(model in unique(aux$Model)){
                if(nrow(aux[aux$Model==model,])==3 && identical(aux[aux$Model==model & aux$Covariate=="cAGE","Parameter"],"varphi[S]") && identical(aux[aux$Model==model & aux$Covariate=="RACE","Parameter"],"varphi[L]") && identical(aux[aux$Model==model & aux$Covariate=="SEX","Parameter"],"delta[Ab]")){
                  NoFNModel = NoFNModel +1
                }
              }
              aux = CovariateModelSelection[CovariateModelSelection$ProjectNumber==paste0(sim," covariates") & CovariateModelSelection$TypeOfSim==type & CovariateModelSelection$Method==meth,]

              TrueModel = 0
              for(model in unique(aux$Model)){
                if(setequal(aux[aux$Model==model,"Covariate"],c("cAGE","SEX","RACE")) && identical(aux[aux$Model==model & aux$Covariate=="cAGE","Parameter"],"varphi[S]") && identical(aux[aux$Model==model & aux$Covariate=="RACE","Parameter"],"varphi[L]") && identical(aux[aux$Model==model & aux$Covariate=="SEX","Parameter"],"delta[Ab]")){
                  TrueModel = TrueModel + 1
                }
              }
              resultModel <- rbind(resultModel,data.frame(ProjectNumber = latex2exp::TeX(paste0(sim," covariates"),output="character"),
                                                          TypeOfSim = type,
                                                          Method = meth,
                                                          NoFNModel = NoFNModel/TotalNumberofModel[[paste0(meth,sim,type)]],
                                                          TrueModel = TrueModel/TotalNumberofModel[[paste0(meth,sim,type)]]))
            }
          }
        }
      }
      CovariateModelSelection[CovariateModelSelection$Method=="reg","Method"] <- "regression"
      resultCovariatePar[resultCovariatePar$Method=="reg","Method"] <- "regression"
      resultCovariate[resultCovariate$Method=="reg","Method"] <- "regression"
      resultModel[resultModel$Method=="reg","Method"] <- "regression"


      resultCovariate[resultCovariate$Method=="CoVSS","Method"] <- "ClustOfVar"
      resultCovariatePar[resultCovariatePar$Method=="CoVSS","Method"] <- "ClustOfVar"
      CovariateModelSelection[CovariateModelSelection$Method=="CoVSS","Method"] <- "ClustOfVar"
      resultModel[resultModel$Method=="CoVSS","Method"] <- "ClustOfVar"

      save(resultCovariate,CovariateModelSelection,orderList,resultCovariatePar,resultModel,file=paste0("Save/BuildParResults_",Rv,".RData"))
    }
  }
}
