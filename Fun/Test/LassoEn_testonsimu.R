library(slurmR)

# Arguments of batch
i <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

source("Fun/Rsmlx/elasticnetSelection.R")
source("Fun/Rsmlx/lassoSelection.R") # classique : intercept=TRUE, all in one, scale Y 
load("Fun/Test/DataSet.RData")

# H0 et H1
H1 =  paste0("(",which(Beta!=0,arr.ind = TRUE)[,"row"],";",which(Beta!=0,arr.ind = TRUE)[,"col"],")") #H1 correspond au coefficient non nul, i.e. ceux qui doivent être sélectionné donc les vrais positives le reste sont des vrais négatif CEUX QU'ON VEUT SELECTIONNER
H0 = setdiff(paste0("(",rep(1:n,each=m),";",1:m,")"),H1) # CEUX QU'ON VEUT PAS SELECTIONNER
cat(paste0("#H1 =",length(H1),".\n")) #21
cat(paste0("#H0 =",length(H0),".\n")) #379

TP <- function(sel,H0,H1){length(unique(intersect(sel,H1)))}
TN <- function(sel,H0,H1){length(unique(setdiff(H0,sel)))}
FP <- function(sel,H0,H1){length(unique(intersect(sel,H0)))}
FN <- function(sel,H0,H1){length(unique(setdiff(H1,sel)))}

FDR <- function(sel,H0,H1){round(FP(sel,H0,H1)/(max(FP(sel,H0,H1)+TP(sel,H0,H1),1)),digits=3)}
FNR <- function(sel,H0,H1){round(FN(sel,H0,H1)/(max(FN(sel,H0,H1)+TN(sel,H0,H1),1)),digits=3)}

Res = data.frame()


for(corrInX in c(1,2,3)){  # 1 pas de corr, 2 corr
  for(corrInY in c(1,2)){
    eval(parse(text=paste0("X <- X_",corrInX)))
    eval(parse(text=paste0("Omega <- Omega_",corrInY)))
    
    res = lapply(1:m,FUN=function(x){lassoSelection(Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,x],X,Omega[x,x])})
    aux = Reduce(rbind,lapply(1:m,FUN=function(x){res[[x]]$selection}))
    sel <- paste0("(",which(aux!=0,arr.ind = TRUE)[,"col"],";",which(aux!=0,arr.ind = TRUE)[,"row"],")")
    
    Res <- rbind(Res,data.frame(Ind = i,
                                alpha=NA,
                                lambda=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"lambda")],digits=3)}),collapse=", "),
                                thresholds = paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"thresholds")],digits=3)}),collapse=", "),
                                corrInX = corrInX,
                                corrInY=corrInY,
                                Method = "LassoCrit",
                                FP = FP(sel,H0,H1),
                                FN = FN(sel,H0,H1),
                                FDR=FDR(sel,H0,H1),
                                FNR=FNR(sel,H0,H1)))
    
    res = lapply(1:m,FUN=function(x){lassoSelection(Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,x],X,Omega[x,x],thresholdsSS = 0.95)})
    aux = Reduce(rbind,lapply(1:m,FUN=function(x){res[[x]]$selection}))
    sel <- paste0("(",which(aux!=0,arr.ind = TRUE)[,"col"],";",which(aux!=0,arr.ind = TRUE)[,"row"],")")
    
    Res <- rbind(Res,data.frame(Ind = i,
                                alpha=NA,
                                lambda=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"lambda")],digits=3)}),collapse=", "),
                                thresholds = paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"thresholds")],digits=3)}),collapse=", "),
                                corrInX = corrInX,
                                corrInY=corrInY,
                                Method = "Lasso",
                                FP = FP(sel,H0,H1),
                                FN = FN(sel,H0,H1),
                                FDR=FDR(sel,H0,H1),
                                FNR=FNR(sel,H0,H1)))
    
    res = lapply(1:m,FUN=function(x){elasticnetSelection(Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,x],X,Omega[x,x])})
    aux = Reduce(rbind,lapply(1:m,FUN=function(x){res[[x]]$selection}))
    sel <- paste0("(",which(aux!=0,arr.ind = TRUE)[,"col"],";",which(aux!=0,arr.ind = TRUE)[,"row"],")")
    
    Res <- rbind(Res,data.frame(Ind = i,
                                alpha=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"alpha")],digits=3)}),collapse=", "),
                                lambda=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"lambda")],digits=3)}),collapse=", "),
                                thresholds = paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"thresholds")],digits=3)}),collapse=", "),
                                corrInX = corrInX,
                                corrInY=corrInY,
                                Method = "ElasticNetCrit",
                                FP = FP(sel,H0,H1),
                                FN = FN(sel,H0,H1),
                                FDR=FDR(sel,H0,H1),
                                FNR=FNR(sel,H0,H1)))
    
    res = lapply(1:m,FUN=function(x){elasticnet(Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,x],X,Omega[x,x],thresholdsSS = 0.95)})
    aux = Reduce(rbind,lapply(1:m,FUN=function(x){res[[x]]$selection}))
    sel <- paste0("(",which(aux!=0,arr.ind = TRUE)[,"col"],";",which(aux!=0,arr.ind = TRUE)[,"row"],")")
    
    Res <- rbind(Res,data.frame(Ind = i,
                                alpha=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"alpha")],digits=3)}),collapse=", "),
                                lambda=paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"lambda")],digits=3)}),collapse=", "),
                                thresholds = paste0(sapply(1:m,FUN=function(x){round(res[[x]]$param[stringr:::str_detect(names(res[[x]]$param),"thresholds")],digits=3)}),collapse=", "),
                                corrInX = corrInX,
                                corrInY=corrInY,
                                Method = "ElasticNet",
                                FP = FP(sel,H0,H1),
                                FN = FN(sel,H0,H1),
                                FDR=FDR(sel,H0,H1),
                                FNR=FNR(sel,H0,H1)))
    
    res = data.frame()
    for(l in 1:m){
      y <- Y.list[[paste0("Y",corrInX,corrInY)]][[i]][,l]
      l_data = cbind.data.frame(y,X)
      rownames(l_data) <- 1:nrow(l_data)
      f.sature <- as.formula("y ~ .")
      model.sature = lm(f.sature, l_data)
      f.cst <- as.formula("y~1")
      model.cst = lm(f.cst, l_data)
      
      llk = lm.cst = stepAIC(model.cst,  scope = list(upper = model.sature,lower = model.cst),trace = FALSE)
      
      res = rbind(res,as.logical(tabulate(as.numeric(stringr::str_remove_all(setdiff(names(llk$coefficients),"(Intercept)"),"Gen")),n)))
    }
    sel <- paste0("(",which(res!=0,arr.ind = TRUE)[,"col"],";",which(res!=0,arr.ind = TRUE)[,"row"],")")
  
    Res <- rbind(Res,data.frame(Ind = i,
                                alpha=NA,
                                lambda=NA,
                                thresholds = NA,
                                corrInX = corrInX,
                                corrInY=corrInY,
                                Method = "stepAIC",
                                FP = FP(sel,H0,H1),
                                FN = FN(sel,H0,H1),
                                FDR=FDR(sel,H0,H1),
                                FNR=FNR(sel,H0,H1)))
    
  }
}

save(Res,file=paste0("Fun/Test/Results/Res_",i,".RData"))