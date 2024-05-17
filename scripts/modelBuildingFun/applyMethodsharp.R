applyMethodsharp <- function(Y,X,omega,cov0,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             criterion="BIC",ncrit=20,covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1,
                             iter=1){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètre
  
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
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){ 
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1
  }else{
    VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,exclude=exclude,nfolds=nfolds,alpha=alpha,K=nSS,n_cores=n_cores)
    selection = sharp::SelectedVariables(VariableSelection.outputs)
  }
  
  plot=CalibrationPlot(VariableSelection.outputs) + ggplot2::ggtitle(paste0("Calibration Plot at iteration ",iter," for parameters ",p.name))
  
  
  if(all(!as.logical(selection))){
    newcriterion = critFUN(lm(Ywh ~ NULL))
  }else{
    Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
    newcriterion = critFUN(lm(Ywh~Xkeep))
  }
  if(oldCriterion == -Inf){
    oldCriterion = Inf
  }
  
  if(newcriterion >= oldCriterion){
    to.cat <- c(to.cat,paste0("        No model improving the criterion as been find, the previous covariate model is kept."))
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,paste0("\n              > parameter values : ",
        paste0(c("lambda","thresholds"),"=",c("","0."),
               sapply(sharp::ArgmaxId(VariableSelection.outputs),function(x){
                 10**floor(log10(x))*round(x/(10**floor(log10(x))),digits=2)
               }),collapse=", ")))
  }
  
  model.list = modelFromSelection(Y,X,selection)
  
  
  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat,plot=list(plot=plot,p.name=p.name)))
}

CalibrationPlot <- function(VariableSelection.outputs){
  pi_list = VariableSelection.outputs$params$pi_list
  lambda_list = VariableSelection.outputs$Lambda
  
  
  Score = VariableSelection.outputs$S_2d
  df = data.frame()
  for(i in 1:ncol(Score)){
    df <- rbind(df,data.frame(lambda = signif(lambda_list,digits=2),pi = pi_list[i],Score=Score[,i]))
  }
  
  plot = ggplot2::ggplot(df,ggplot2::aes(x=as.factor(lambda),y=pi,fill=Score)) +
    ggplot2::geom_tile() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggplot2::scale_fill_gradientn(colours = c("ivory", "navajowhite", "tomato","darkred"),na.value="white") + 
    ggplot2::xlab(latex2exp::TeX("$\\lambda$"))+ggplot2::ylab(latex2exp::TeX("$\\pi$")) +
    ggplot2::theme(axis.title=ggplot2::element_text(size=14,face="bold")) +
    ggplot2::scale_y_continuous(breaks=c(seq(0,0.5,0.05),seq(0.6,1,0.1)))
  
  return(plot)
}
