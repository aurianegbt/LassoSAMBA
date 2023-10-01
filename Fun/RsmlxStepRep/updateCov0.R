## Cette fonction prends en entrée une matricede covariables et d'individus X, Y tel que : 
##     -  Y est un vecteur donc pour UN paramètre avec en ligne les différents individus 
##     - Y est une matrice avec sur chaque colonne le paramètre 
##     
## Sigma est la matrice de covariance des paramètres considérés 
## X est la matrice des covariables 
##     
## A ce niveau là, rien n'est BLANCHI ni STANDARDIZE, X et Y sont des matrices quantitatives. 
## cov0 correspond au colonne de X a exclure de la sélection (doit être du même type)

updateCov0 <- function(Y, eta.list, x, p.max=1, covFix = NULL,
                       pen.coef = NULL, pw.list = NULL, cov0.list = NULL){
  param = names(cov0.list)
  n.param=length(param)
  t.param = setdiff(colnames(Y),c("id","rep"))
  for(j in 1:n.param){
    cov0.list[[param[j]]] <- upd(param[j],Y[,t.param[j],drop=FALSE],eta.list[[param[j]]],x,p.max,covFix,pen.coef,pw.list[[param[j]]],cov0.list[[param[j]]])
  }
  return(cov0.list)
}

upd <- function(ny,y,eta,x,p.max,covFix,pen.coef,pw,cov0){
  N <- length(unique(x$id))
  nrep <- nrow(y)/N
  if (p.max <= 1) {
    nx <- setdiff(names(x), c("id", "rep"))
    yc <- rowMeans(matrix(y[[1]], nrow = N))
    xc <- x
    lm0 <- lm(yc ~ 1)
    nxc <- setdiff(nx, cov0)
    pjc <- NULL
    for (nc in nxc) {
      lmc <- lm(yc ~ xc[[nc]])
      pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
      pjc <- c(pjc, pc)
    }
    pjc <- Rsmlx:::p.weight(pjc, pw[nxc], pen.coef)
    names(pjc) <- nxc
    if (!is.null(eta)) {
      etac <- rowMeans(matrix(eta[[1]], nrow = N))
      lm0 <- lm(etac ~ 1)
      pjec <- NULL
      for (nc in nxc) {
        lmc <- lm(etac ~ xc[[nc]])
        pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
        pjec <- c(pjec, pc)
      }
      pjec <- Rsmlx:::p.weight(pjec, pw[nxc], pen.coef)
      names(pjec) <- nxc
      pjc <- pmin(pjc, pjec, na.rm = T)
      
    }
    pjc[names(which(Rsmlx:::mlx.getIndividualParameterModel()$covariateModel[[ny]]))] <- 0
    list.c <- which(pjc > p.max)
    cov0 <- c(cov0, nxc[list.c])
  }
  return(cov0)
}
