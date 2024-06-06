sharpCalibration <- function(Y,X,omega,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             p.name=NULL,
                             n_cores=NULL,
                             iter = 1){
  
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
  
  VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,nfolds=nfolds,alpha=alpha,K=nSS,n_cores=n_cores)
  
  pi_list = VariableSelection.outputs$params$pi_list
  lambda_list = VariableSelection.outputs$Lambda
  
  
  Score = VariableSelection.outputs$S_2d
  df = data.frame()
  for(i in 1:ncol(Score)){
    df <- rbind(df,data.frame(lambda = signif(lambda_list,digits=2),pi = pi_list[i],Score=Score[,i]))
  }
  
  
  CIl = range(df[df$Score>=quantile(df$Score,0.90,na.rm=T),"lambda"],na.rm=T)
  CIpi = quantile(df[df$Score>=quantile(df$Score,0.90,na.rm=T),"pi"],0.90,na.rm=T) 
  
  plot = ggplot(df,aes(x=as.factor(lambda),y=pi,fill=Score)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_gradientn(colours = c("ivory", "navajowhite", "tomato","darkred"),na.value="white") + 
    xlab(latex2exp::TeX("$\\lambda$"))+ylab(latex2exp::TeX("$t_{SS}$")) +
    theme(axis.title=element_text(size=14,face="bold")) +
    theme(axis.text.x=element_text(size=8)) +
    scale_y_continuous(breaks=c(seq(0,0.5,0.05),seq(0.6,1,0.1)))  + 
    ggtitle(latex2exp::TeX(r"(Heatmap of the stability score S$(t_{SS},\lambda)$)")) + 
    theme(plot.subtitle = element_text(size=10)) + 
    theme(legend.position = "bottom") + 
    theme(legend.text = element_text(face="plain"))  +
    theme(legend.title = element_text(size=10))
  
  return(list(lambda.grid = seq(min(CIl),max(CIl),length.out =100),thresholdsSS=unname(CIpi),p.name=p.name,plot=plot))
}