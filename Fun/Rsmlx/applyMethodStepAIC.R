applyMethodStepAIC <- function(Y,X,
                             omega,cov0,cov1,covariate.model=NULL,n.full=10,
                             direction="both",steps=1000,nb.model=1,
                             pweight=NULL,pen.coef,p.name=NULL){
  # X et Y on juste la bonne "forme" mais pas scale encore (ça c'est fait dans la sélection en elle même !!!! )
  # Le but de cette fonction est de construit le modèle linéaire pour chaque paramètr
  cov.names = colnames(X)
  tparam.names = colnames(Y)

  resSelection = StepAICSelection(Y,X,omega,cov0,cov1,covariate.model,n.full,direction,steps,nb.model,pweight,pen.coef,p.name)

  return(list(model=resSelection$model,res=resSelection$selection,cov0=cov0,p.name=p.name))
}
