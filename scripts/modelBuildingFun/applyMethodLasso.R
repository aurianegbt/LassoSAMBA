applyMethodLasso <- function(Y,X,omega,cov0,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             criterion="BIC",
                             covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1,
                             iter=1,
                             FDR_thr=0.10){
  to.cat = c()
  
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
  if(!is.matrix(X)){
    Xaux <- as.matrix(X)
  }else{Xaux=X}
  
  Xsc <- scale(apply(Xaux,2,FUN=as.numeric))
  
  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  Ywh <- Yaux %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)
  colnames(Xwh) <- cov.names

  if(!is.null(covariate.model)){
    savedSelection = setNames(as.numeric(covariate.model),names(covariate.model))
    prevSelection = covariate.model
    
    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(prevSelection)[which(prevSelection)]]
      oldCriterion = critFUN(lm(Ywh ~ Xkeep))
    }
    to.cat <- c(to.cat,paste0("\n Lasso selection, calibrated using sharp method, improving the ",criterion," criterion for ",p.name," :\n "))
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
      newcriterion = critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
      newcriterion = critFUN(lm(Ywh~Xkeep))
    }
    
    to.cat.here = ""
  }else{
    VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,exclude=exclude,nfolds=nfolds,pi_list=seq(0.50,0.99,0.01),alpha=alpha,K=nSS,n_cores=n_cores,FDP_thr = FDR_thr)
    
    pi_list = VariableSelection.outputs$params$pi_list
    lambda_list = VariableSelection.outputs$Lambda
    
    
    Score = VariableSelection.outputs$S_2d
    argmax_id = which(!is.na(Score),arr.ind=T)
    if(nrow(argmax_id)!=0){
      resSharp = lapply(split(argmax_id,1:nrow(argmax_id)),FUN=function(arg_id){
        selection = sharp::SelectedVariables(VariableSelection.outputs,argmax_id = arg_id)
        if(all(!as.logical(selection))){
          newcriterion = critFUN(lm(Ywh ~ NULL))
        }else{
          Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
          newcriterion = critFUN(lm(Ywh~Xkeep))
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
      
      plot1 =
        ggplot2::ggplot(df,ggplot2::aes(x=as.factor(lambda),y=pi,fill=Score)) +
        ggplot2::geom_tile() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        ggplot2::scale_fill_gradientn(colors = c("ivory", "navajowhite", "tomato","darkred","black"),na.value="white") + 
        ggplot2::xlab(latex2exp::TeX("$\\lambda$"))+
        ggplot2::ylab(latex2exp::TeX("$t_{SS}$")) +
        ggplot2::theme(axis.title=ggplot2::element_text(size=14,face="bold")) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(size=8)) +
        ggplot2::scale_y_continuous(breaks=c(seq(0,0.5,0.05),seq(0.6,1,0.1)))  + 
        ggplot2::ggtitle(latex2exp::TeX(r"(Heatmap of the constrained stability score S$(t_{SS},\lambda)$ )")) + 
        ggplot2::theme(plot.subtitle = ggplot2::element_text(size=10)) + 
        ggplot2::theme(legend.position = "bottom") + 
        ggplot2::theme(legend.text = ggplot2::element_text(face="plain"))  +
        ggplot2::theme(legend.title = ggplot2::element_text(size=10))
      
      to.cat.here =  paste0("\n              > parameter values : ",
                            paste0(c("lambda","thresholds"),"=",c(signif(lambda_list[argmax_id[indMax,1]],3),signif(pi_list[argmax_id[indMax,2]],2)),collapse=", "))
    }else{
      selection = savedSelection
      newcriterion = oldCriterion
      
      plot1=NULL
    }
  }
  
  if(newcriterion >= oldCriterion){
    to.cat <- c(to.cat,paste0("        No model improving the criterion as been find, the previous covariate model is kept."))
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,to.cat.here)
  }
  
  model.list = modelFromSelection(Y,X,selection)
  
  
  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat,plot=list(plot1=plot1)))
}