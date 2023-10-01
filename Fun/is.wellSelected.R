is.wellSelected <- function(Model=NULL){
  if(is.null(Model))
    Model <- Rsmlx:::mlx.getIndividualParameterModel()
  
  flag = TRUE
  # Statistical Model :
  if(flag){
    flag <- all(Model$distribution=="logNormal")
    flag <- all(Model$variability[[1]] == c(rep(FALSE,2),rep(TRUE,3)))
  }
  # Structural Model :
  if(flag){
    trueModel = list(delta_S=character(0),delta_L=character(0),phi_S="cAGE",phi_L="RACE",delta_AB="SEX")
    CovariateModel = list()
    for (parm in names(Model$covariateModel)){
      aux = Model$covariateModel[[parm]]
      CovariateModel <- append(CovariateModel,list(names(aux[aux])))
      names(CovariateModel)[length(CovariateModel)] <- parm
    }
    flag <- identical(CovariateModel,trueModel)
  }
  return(flag)
}