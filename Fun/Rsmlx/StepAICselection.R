StepAICSelection <- function(Y,X,
                           omega,cov0,cov1,covariate.model=NULL,n.full=10,
                           direction="both",steps=1000,nb.model=1,pweight=NUL,pen.coef=NULL,p.name=NULL){

  cov.names = colnames(X)

  l_data = cbind.data.frame(y=Y[[1]],X)
  
  if(length(setdiff(cov.names,cov0))<=n.full){
    direction = "full"
  }
  
  if(direction!="full"){
    llk = tryCatch({
      f.sature <- "y ~ ."
      if (length(cov0) > 0) 
        f.sature <- paste0(f.sature, "-", paste(cov0, 
                                                collapse = "-"))
      f.sature <- as.formula(f.sature)
      model.sature = lm(f.sature, l_data)
      f.cst <- "y ~ 1"
      if (length(cov1) > 0) 
        f.cst <- paste0(f.cst, "+", paste(cov1, collapse = "+"))
      f.cst <- as.formula(f.cst)
      model.cst = lm(f.cst, l_data)
      if (direction == "backward") {
        lm.sat = MASS::stepAIC(model.sature, direction = "backward", 
                             trace = FALSE, k = nrep * pen.coef, scope = list(upper = model.sature, 
                                                                              lower = model.cst), steps = steps, weight = pweight)
      } else {
        c.cur <- names(which(Rsmlx:::mlx.getIndividualParameterModel()$covariateModel[[p.name]]))
        f.cur <- "y ~ 1"
        if (length(c.cur) > 0) 
          f.cur <- paste0(f.cur, "+", paste(c.cur, collapse = "+"))
        f.cur <- as.formula(f.cur)
        model.cur = lm(f.cur, l_data)
        lm.cst = MASS::stepAIC(model.cur, direction = direction, 
                         trace = FALSE, k = nrep * pen.coef, scope = list(upper = model.sature, 
                                                                          lower = model.cst), steps = steps)
      }
    }, error = function(e) {
      print("Error in stepAIC")
      return(-Inf)
    })
    Gnames = colnames(X)
    G <- data.frame(matrix(0, ncol = length(Gnames), nrow = 1))
    colnames(G) <- Gnames
    usedcovariates = names(llk$model)[-1]
    G[1, usedcovariates] <- 1
    ll = logLik(llk)/nrep
    df = length(coef(llk)) - 1
    crit = -2 * ll + pen.coef * sum(pweight[usedcovariates])
    res <- data.frame(ll = round(ll, digits = 3), df = df, 
                      criterion = round(crit, digits = 3))
    
    res <- cbind(G == 1, res)
    if (nb.model == 1) 
      res[, c("ll", "df", "criterion")] <- NULL
    else res[2:nb.model, ] <- res[1, ]
    row.names(res) <- 1:nrow(res)
    return(list(model = llk, selection=res))
  }else{
    s <- rep("0:1", length(cov.names))
    names(s) <-cov.names
    if (!is.null(cov0)) 
      s[cov0] <- "0"
    s <- paste(s, collapse = ",")
    s <- paste0("G <- expand.grid(", s, ")")
    eval(parse(text = s))
    names(G) <- cov.names
    if (length(cov0) > 0) {
      i0 <- which(rowSums(G[cov0]) == 0)
      G <- G[i0, ]
    }
    if (length(cov1) > 0) {
      i1 <- which(rowSums(G[cov1] == 1) == length(cov1))
      G <- G[i1, ]
    }
    ng <- nrow(G)
    d <- ncol(G)
    ll <- df <- bic  <- NULL
    yg <- Y[[1]]
    xg <- X
    for (k in 1:ng) {
      xk <- data.frame(y = yg)
      Gk <- G[k, , drop = FALSE]
      pwk <- pweight[names(Gk)]
      j1 <- which(Gk == 1)
      if (length(j1) > 0) 
        xk[names(X)[j1]] <- xg[j1]
      llk = tryCatch({
        lmk <- lm(y ~ ., data = xk)
        logLik(lmk)[1]/nrep
      }, error = function(e) {
        return(-Inf)
      })
      dfk <- sum(Gk > 0)
      bick <- -2 * llk + pen.coef * sum(Gk * pwk)
      ll <- c(ll, llk)
      df <- c(df, dfk)
      bic <- c(bic, bick)
    }
    bic <- round(bic, digits = 3)
    i0 <- rep(1, ng)
    mG <- ncol(G)
    for (k in seq_len(ng - 1)) {
      if (i0[k] == 1) {
        ik <- which(bic[(k):ng] == bic[k]) + k - 1
        sk <- .rowSums(G[ik, ] == 2, n = length(ik), m = mG)
        ik0 <- ik[which(sk == 0)]
        if (length(ik0) == 0) 
          ik0 <- ik[order(sk)[1]]
        i0[ik] <- 0
        i0[ik0] <- 1
        i0[k] <- 1
      }
    }
    res <- data.frame(ll = round(ll, digits = 3), df = df, criterion = bic)
    res <- res[i0 == 1, ]
    G <- G[i0 == 1, , drop = FALSE]
    bic <- bic[i0 == 1]
    eval(parse(text = paste0(names(Y), " <- Y[[1]]")))
    obic <- order(bic)
    k.min <- obic[1]
    Gkmin <- G[k.min, ]
    j1 <- which(Gkmin == 1)
    if (length(j1) > 0) {
      for (k in (1:length(j1))) eval(parse(text = paste0(names(X)[j1[k]], 
                                                         " <- rep(X[[j1[k]]], nrep)")))
    }
    list.x <- c("1", names(X)[j1])
    form1 <- paste0(names(Y), "~", paste(list.x, collapse = "+"))
    eval(parse(text = paste0("lm.min <- lm(", form1, ")")))
    lm.min$covsel = Gkmin
    nb.model0 <- min(nb.model, length(bic))
    res <- res[obic[1:nb.model0], ]
    G <- G[obic[1:nb.model], , drop = FALSE]
    if (nb.model > nb.model0) 
      res[(nb.model0 + 1):nb.model, c("ll", "df", "criterion")] <- NA
    res <- cbind(G == 1, res)
    if (nb.model == 1) 
      res[, c("ll", "df", "criterion")] <- NULL
    row.names(res) <- 1:nrow(res)
    return(list(model = lm.min, selection=res))
  }
}
