library(simstudy)

# Dataset generation
set.seed(1710)
n = 200
m = 4
N = 100

load("Fun/Test/correlationMatrix.RData")

param_pop = rnorm(m)
Y_0 = matrix(rep(param_pop,N),byrow = T,nrow=N)

Beta = matrix(rbinom(n*m,1,0.05),nrow=n)*rnorm(n*m)

Sigma_1 = diag(rexp(n,1.5))
Sigma_2 = genCorMat(n)
ind=sample(1:ncol(corMatrix),size = n,replace = F)
Sigma_3 = corMatrix[ind,ind]

mu = rnorm(n)

X_1 = MASS::mvrnorm(N,mu,Sigma_1)
X_2 = MASS::mvrnorm(N,mu,Sigma_2)
X_3 = MASS::mvrnorm(N,mu,Sigma_3)
colnames(X_1) <- colnames(X_2) <- colnames(X_3) <- paste0("Gen",1:n)

Omega_1 = diag(rexp(m,2.5))
Omega_2 = genCorMat(m)

Y.list = list(Y11=list(),Y12=list(),Y21=list(),Y22=list(),Y31=list(),Y32=list())

for(i in 1:50){
  eta_1 = Reduce(rbind,lapply(1:N,function(i){MASS::mvrnorm(1,rep(0,m),Omega_1)}))
  eta_2 = Reduce(rbind,lapply(1:N,function(i){MASS::mvrnorm(1,rep(0,m),Omega_2)}))
  
  Y_11 = Y_0 + X_1 %*% Beta + eta_1
  Y.list$Y11 <- append(Y.list$Y11,list(Y_11))
  
  
  Y_12 = Y_0 + X_1 %*% Beta + eta_2
  Y.list$Y12 <- append(Y.list$Y12,list(Y_12))
  
  Y_21 = Y_0 + X_2 %*% Beta + eta_1
  Y.list$Y21 <- append(Y.list$Y21,list(Y_21))
  
  
  Y_22 = Y_0 + X_2 %*% Beta + eta_2
  Y.list$Y22 <- append(Y.list$Y22,list(Y_22))
  
  Y_31 = Y_0 + X_3 %*% Beta + eta_1
  Y.list$Y31 <- append(Y.list$Y31,list(Y_31))
  
  
  Y_32 = Y_0 + X_3 %*% Beta + eta_2
  Y.list$Y32 <- append(Y.list$Y32,list(Y_32))
  
}

save(n,m,N,Y_0,Y.list,Beta,Omega_1,Omega_2,X_1,X_2,X_3,file="Fun/Test/DataSet.RData")
