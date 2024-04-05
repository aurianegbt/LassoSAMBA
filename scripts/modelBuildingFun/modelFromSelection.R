modelFromSelection <- function(Y,X,selection=NULL){
  if(!is.data.frame(Y)){
    Y <- as.data.frame(Y)
  }
  if(!is.data.frame(X)){
    X <- as.data.frame(X)
  }

  m = ncol(Y)
  lm.list=list()
  for(i in 1:m){
    data = cbind(y = matrix(Y[,i],ncol=1),X)
    formula = paste0("y ~ 1 ")
    if(m==1){list.c = names(selection)[selection==1]}else{list.c = colnames(selection)[selection[i,]==1]}
    if(length(list.c)>0){
      formula = paste0(formula," + ",paste0(list.c,collapse=" + "))
    }
    eval(parse(text=paste0("lm.sel <- lm(",formula,",data=data)")))
    lm.list <- append(lm.list,list(lm.sel))
  }

  if(m==1){
    return(lm.sel)
  }else{
    return(lm.list)
  }
}
