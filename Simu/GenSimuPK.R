library(dplyr)
library(foreach)
set.seed(41520810)
dir <- function(d){if(!dir.exists(d)){dir.create(d)}}
dir("Files/FilesPK")
dir("Files/FilesPK/200cov")
dir("Files/FilesPK/200corcov")
dir("Files/FilesPK/200cov/simulation");  dir("Files/FilesPK/200cov/covTable")
dir("Files/FilesPK/200corcov/simulation");  dir("Files/FilesPK/200corcov/covTable")
# load("Simu/DataTransCoding.RData")
# 
# aux = (dataTransCoding %>% filter(visit=="M6"))[,c(7,9:18662)]
# genesKept = names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[2:201][sample(1:200)]
# 
# aux <- aux[,genesKept]
# colnames(aux) <- c(1:ncol(aux))
# 
# genCorMat <- cor(aux,method="spearman")
# epsilon <- 1e-10
# genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))
# 
# save(genCorMat,aux,file="Simu/corrPK.RData")
load("Simu/corrPK.RData")

N=100

ka_pop = 0.1
C0_pop = 2000


mu = rep(0,ncol(aux))
sd = rep(1,ncol(aux))

cov = 1:5
Beta = rnorm(n=length(cov),mean=0,sd=1)

omega_ka = 1
omega_C0 = 1
sigma = 1

t_ij <-c(seq(0,30,1),50,100)

tocheck = list(numeric(),numeric())
tocheck_corr = numeric()

for(r in 1:100){
  covariates = mvtnorm::rmvnorm(n=N,mean = mu,sigma=diag(sd))
  covariates_corr = mvtnorm::rmvnorm(n=N,mean = mu,sigma = genCorMat)
  
  colnames(covariates) <- colnames(covariates_corr) <- paste0("Gen",1:200)
  
  ka_i = exp(log(ka_pop)+covariates[,cov]%*%Beta+rnorm(N,0,omega_ka))
  C0_i = exp(log(C0_pop)+rnorm(N,0,omega_C0))
  
  tocheck[[1]] <- c(tocheck[[1]],ka_i)
  tocheck[[2]] <- c(tocheck[[2]],C0_i)
  
  ka_i_corr = exp(log(ka_pop)+covariates_corr[,cov]%*%Beta+rnorm(N,0,omega_ka))
  
  tocheck_corr <- c(tocheck_corr,ka_i_corr)

  sim = data.frame()
  sim_corr = data.frame()

  for(i in 1:N){
    sim <- rbind(sim,
                 data.frame(id=i,
                            time=t_ij,
                            yC=C0_i[i]*exp(-ka_i[i]*t_ij) + (0.1*C0_i[i]*exp(-ka_i[i]*t_ij) + 20)*rnorm(length(t_ij),0,sigma)))

    sim_corr <- rbind(sim,
                      data.frame(id=i,
                                 time=t_ij,
                                 yC=C0_i[i]*exp(-ka_i_corr[i]*t_ij) + (0.1*C0_i[i]*exp(-ka_i_corr[i]*t_ij) + 20)*rnorm(length(t_ij),0,sigma)))
  }
  sim <- merge(sim,cbind(id=1:N,covariates),by="id")
  sim_corr <- merge(sim_corr,cbind(id=1:N,covariates_corr),by="id")
  
  
  
  write.csv(sim,file=paste0("Files/FilesPK/200cov/simulation/simulation_",r,".txt"),quote = F,row.names = F)
  write.csv(sim_corr,file=paste0("Files/FilesPK/200corcov/simulation/simulation_",r,".txt"),quote = F,row.names = F)
  if(r==1){
   covTable <- cbind(id=1:N,covariates) 
   corcovTable <- cbind(id=1:N,covariates_corr) 
   write.csv(covTable,file=paste0("Files/FilesPK/200cov/covTable/covTable_",r,".txt"),quote = F,row.names = F)
   write.csv(corcovTable,file=paste0("Files/FilesPK/200corcov/covTable/covTable_",r,".txt"),quote = F,row.names = F)
  }
}

headerTypes <- c("id","time","observation",rep("contcov",200))
save(headerTypes,file="Files/FilesPK/200cov/headerTypes.RData")
save(headerTypes,file="Files/FilesPK/200corcov/headerTypes.RData")
