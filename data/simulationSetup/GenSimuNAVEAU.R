dir <- function(path){if(!dir.exists(path)){dir.create(path)}}

## Load the library required 
suppressMessages({
  library(dplyr)
  library(lixoftConnectors)
  library(simstudy)
  library(data.table)
  library(ggcorrplot)
  library(ggpubr,quietly=TRUE)
  
  initializeLixoftConnectors("simulx")
})

library(lixoftConnectors)
library(RsSimulx)
library(Rsmlx)
library(dplyr)
data.dir<- "data/simulationFiles/FilesNAVEAU"
dir(data.dir)
data.dir<- "data/simulationFiles/FilesNAVEAU/simulation"
dir(data.dir)
data.dir.cov <- "data/simulationFiles/FilesNAVEAU/covTable"
dir(data.dir.cov)
project_ini.dir <- "project_ini"
project_star.dir <- "project_star"
file.results <- "simulations.RData"

#Simulation settings
M <- 100 #number replicates
seed0 <- 12345

N <- 200 #Number of individuals
ncov <- 500 #Number of covariates

## Generate 200 correlated covariates and then create 100 replicates
t.y <- seq(150, 3000,length=10) #observation times

# Model specification
lib <- 'naveau_simulation.txt'
parm <- c(phi_pop=1200,
            psi1=200,
            psi2=300,
            omega_phi=sqrt(200),
            beta_C1=100,
            beta_C2=50,
            beta_C3=20,
            a=sqrt(30))

model <- inlineModel("
[COVARIATE]
input = {C1, C2, C3}

[INDIVIDUAL]
input = {phi_pop,omega_phi,C1,C2,C3,beta_C1,beta_C2,beta_C3}

DEFINITION:
phi = {distribution=normal, typical=phi_pop, covariate={C1,C2,C3}, coefficient={beta_C1,beta_C2,beta_C3}, sd=omega_phi}

[LONGITUDINAL]
input = {phi,psi1,psi2,a}

EQUATION:
Y = psi1 /(1+exp(-(t-phi)/psi2))

DEFINITION:
yY = {distribution=normal, prediction=Y, sd=a}
")

p.pop <- parm[1:3]
p.name <- c("phi", "psi1", "psi2")
structural.model <- lib
model <- model

nparam <- length(p.name)
cov.name <- paste0("C",1:ncov)

data.files <- paste0("simulation_",1:M,".txt")
project_ini.files <- paste0("project_ini",c(paste0("0",1:9),10:M),".mlxtran")
project_star.files <- paste0("project_star",c(paste0("0",1:9),10:M),".mlxtran")

seed <- seed0 + 1:M

# ------------------------------------------------
# Data generation
for (m in 1:M) {
  print(m)
  set.seed(seed[m])
  covariate <- data.frame(id=factor(1:N), mvtnorm::rmvnorm(N,mean=rep(0,ncov)))
  names(covariate)[-1] <- cov.name
  if(m==1){
    write.csv(covariate, file=paste0(data.dir.cov,"/covTable_1.txt"),quote = F,row.names = F)
  }
  
  
  r <- simulx(model     = model, 
              parameter = parm,
              covariate = covariate[c('id','C1', 'C2','C3')],
              output    = list(name="yY", time=t.y),
              settings  = list(seed=seed[m]))
  
  d <- r$yY
  d <- left_join(d, covariate, by="id")
  d[is.na(d)] <- "."
  write.csv(d, file=file.path(data.dir,data.files[m]), quote=F, row.names=F)
}
