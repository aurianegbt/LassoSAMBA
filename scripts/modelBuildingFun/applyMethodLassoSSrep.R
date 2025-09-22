applyMethodLassoSSrep <- function(Y,X,omega,cov0,
                                  Y.means,
                             nfolds=5,
                             alpha=1,
                             criterion="BIC",
                             covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1,
                             iter=1,
                             FDP_thr=0.10){
  
  to.cat <- c()
  
  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC 
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }
  
  cov.names = colnames(X)
  tparam.names = colnames(Y)
  
  if(!is.matrix(Y)){
    Yaux <- as.matrix(Y)
  }else{Yaux=Y}
  if(!is.matrix(Y.means)){
    Yaux.means <- as.matrix(Y.means)
  }else{Yaux.means=Y.means}
  if(!is.matrix(X)){
    Xaux <- as.matrix(X)
  }else{Xaux=X}
  
  Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
  
  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  
  Ywh <- Yaux
  Ywh.means <- Yaux.means %*% rootInvOmega
  Ywh[,2] <- Yaux[,2,drop=FALSE] %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)
  colnames(Xwh) <- cov.names
  
  
  if(!is.null(covariate.model)){
    if(any(!(cov.names %in% names(covariate.model)))){
      savedSelection <- setNames(rep(0,length(cov.names)),cov.names)
      savedSelection[names(covariate.model)] <- as.numeric(covariate.model)
    }else{
      savedSelection = setNames(as.numeric(covariate.model),names(covariate.model))
    }
    prevSelection = covariate.model
    
    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh.means ~ NULL))
    }else{
      Xkeep = Xwh[,names(prevSelection)[which(prevSelection)]]
      oldCriterion = critFUN(lm(Ywh.means ~ Xkeep))
    }
    to.cat <- c(to.cat,paste0("\n Lasso selection, calibrated using sharp method on replicates, improving the ",criterion," criterion for ",p.name," :\n "))
    to.cat <- c(to.cat,paste0("       -> Old Criterion : ",round(oldCriterion,digits=2)),"\n")
  }
  
  if(is.null(cov0)){
    exclude = NULL
  }else{
    exclude = which(cov.names %in% cov0)
  }
  
  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))
    
    to.cat.here = ""
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ 
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1
    if(all(!as.logical(selection))){
      newcriterion = critFUN(lm(Ywh.means ~ NULL))
    }else{
      Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
      newcriterion = critFUN(lm(Ywh.means~Xkeep))
    }
    
    to.cat.here = ""
  }else{
    VariableSelection.outputs = sharp.VariableSelection(Xwh,Ywh,Ywh.means,exclude=exclude,nfolds=nfolds,pi_list=seq(0.50,0.99,0.01),alpha=alpha,K=max(Ywh[,1]),n_cores=n_cores,FDP_thr = FDP_thr)
    
    pi_list = VariableSelection.outputs$params$pi_list
    lambda_list = VariableSelection.outputs$Lambda
    
    
    Score = VariableSelection.outputs$S_2d
    argmax_id = which(!is.na(Score),arr.ind=T)
    if(nrow(argmax_id)!=0){
      resSharp = lapply(split(argmax_id,1:nrow(argmax_id)),FUN=function(arg_id){
        selection = sharp::SelectedVariables(VariableSelection.outputs,argmax_id = arg_id)
        if(all(!as.logical(selection))){
          newcriterion = critFUN(lm(Ywh.means ~ NULL))
        }else{
          Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
          newcriterion = critFUN(lm(Ywh.means~Xkeep))
        }
        
        if(newcriterion==-Inf){
          newcriterion = oldCriterion + 1
        }
        
        return(list(selection=selection,criterion=newcriterion))
      })
      
      indMax = which.min(sapply(resSharp,FUN=function(r){r$criterion}))
      selection = resSharp[[indMax]]$selection
      newcriterion = resSharp[[indMax]]$criterion
      
      df = data.frame()
      for(i in 1:ncol(Score)){
        df <- rbind(df,data.frame(lambda = signif(lambda_list,digits=2),pi = pi_list[i],Score=Score[,i]))
      }
      
      to.cat.here =  paste0("\n              > parameter values : ",
                            paste0(c("lambda","thresholds"),"=",c(signif(lambda_list[argmax_id[indMax,1]],3),signif(pi_list[argmax_id[indMax,2]],2)),collapse=", "))
    }else{
      selection = savedSelection
      newcriterion = oldCriterion
    }
  }
  
  if(newcriterion >= oldCriterion){
    to.cat <- c(to.cat,paste0("        No model improving the criterion as been find, the previous covariate model is kept."))
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,to.cat.here)
  }
  
  model.list = modelFromSelection(Y.means,X,selection)
  
  
  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat))
}