
library(foreach)

doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()))

foreach(i = 1:50) %dopar% {
  library(slurmR)
  library(sharp)
  load("Fun/Test/DataSet.RData")
  
  # H0 et H1
  H1 =  paste0("(",which(Beta!=0,arr.ind = TRUE)[,"row"],";",which(Beta!=0,arr.ind = TRUE)[,"col"],")") #H1 correspond au coefficient non nul, i.e. ceux qui doivent être sélectionné donc les vrais positives le reste sont des vrais négatif CEUX QU'ON VEUT SELECTIONNER
  H0 = setdiff(paste0("(",rep(1:n,each=m),";",1:m,")"),H1) # CEUX QU'ON VEUT PAS SELECTIONNER
  cat(paste0("#H1 =",length(H1),".\n")) 
  cat(paste0("#H0 =",length(H0),".\n")) 
  
  TP <- function(sel,H0,H1){length(unique(intersect(sel,H1)))}
  TN <- function(sel,H0,H1){length(unique(setdiff(H0,sel)))}
  FP <- function(sel,H0,H1){length(unique(intersect(sel,H0)))}
  FN <- function(sel,H0,H1){length(unique(setdiff(H1,sel)))}
  
  FDR <- function(sel,H0,H1){round(FP(sel,H0,H1)/(max(FP(sel,H0,H1)+TP(sel,H0,H1),1)),digits=3)}
  FNR <- function(sel,H0,H1){round(FN(sel,H0,H1)/(max(FN(sel,H0,H1)+TN(sel,H0,H1),1)),digits=3)}
  
  # Res = data.frame()
  Res = get(load(paste0("Fun/Test/Results/Res_",i,".RData")))
  
  for(corrInX in c(1,2,3)){  # 1 pas de corr, 2 corr
    for(corrInY in c(1,2)){
      eval(parse(text=paste0("X <- X_",corrInX)))
      eval(parse(text=paste0("Omega <- Omega_",corrInY)))
      
      res = lapply(1:m,FUN=function(x){sharp::VariableSelection(X,Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,x])})
      lambda = sapply(1:m,FUN=function(x){res[[x]]$Lambda[ArgmaxId(res[[x]])[,1]]})
      thresholds = sapply(1:m,FUN=function(x){res[[x]]$params$pi_list[ArgmaxId(res[[x]])[,2]]})
      
      aux = Reduce(rbind,lapply(1:m,FUN=function(x){sharp::SelectedVariables(res[[x]])}))
      sel <- paste0("(",which(aux!=0,arr.ind = TRUE)[,"col"],";",which(aux!=0,arr.ind = TRUE)[,"row"],")")
      
      Res <- rbind(Res,data.frame(Ind = i,
                                  alpha=NA,
                                  lambda=paste0(sapply(lambda,FUN=function(x){signif(x,4)}),collapse=", "),
                                  thresholds = paste0(thresholds,collapse=", "),
                                  corrInX = corrInX,
                                  corrInY=corrInY,
                                  Method = "SHARP",
                                  FP = FP(sel,H0,H1),
                                  FN = FN(sel,H0,H1),
                                  FDR=FDR(sel,H0,H1),
                                  FNR=FNR(sel,H0,H1)))
    }
  }
  
  save(Res,file=paste0("Fun/Test/Results/Res_",i,".RData"))
}
