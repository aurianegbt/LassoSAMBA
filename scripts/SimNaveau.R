# Empty memory ------------------------------------------------------------
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Load packages, arguments, functions and  data  --------------------------
suppressWarnings({
  suppressMessages({
    library(ggplot2)
    library(nlme)
    library(cowplot)
    library(glmnet)
    library(mvnfast)
    library(doParallel)
    library(dplyr)
  })
})

arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))
sim <- read.csv(paste0("data/simulationFiles/FilesNAVEAU/simulation/simulation_",arr,".txt"))

seed <- arr 
n <- 200 # number of individuals
J <- 10 # number of repetitions
p <- 500 # number of covariates
sigma2 <- 30 # residual variance of y
Gamma2 <- 200 # inter-individual variance of phi_i
mu <- 1200 # intercept
beta <- c(100, 50, 20, rep(0, p - 3)) # covariate fixed effects vector
beta_tilde <- c(mu, beta) # useful reformulation
psi1 <- 200 # parameter of the function g
psi2 <- 300 # parameter of the function g
psi <- c(psi1, psi2)
t <- sim$time
Y <- sim$yY
id <- sim$id
V_tilde <- as.matrix(cbind(rep(1,n),sim[sim$time==150,-c(1:4)]))



betatilde_init <- c(1500, rep(100, 10), rep(1, p - 10))
sigma2_init <- 100
Gamma2_init <- 5000
alpha_init <- 0.5
tau <- 0.98

nu0 <- 0.02 # spike parameter
nu1 <- 12000 # slab parameter
nu_Gamma <- 1 # prior parameter of Gamma2
lb_Gamma <- 1 # prior parameter of Gamma2
nu_sigma <- 1 # prior parameter of sigma2
lb_sigma <- 1 # prior parameter of sigma2
a <- 1 # prior parameter of alpha
b <- p # prior parameter of alpha
sigma2_mu <- 3000^2 # prior parameter of mu

param_init <- list(
  beta_tilde = betatilde_init, alpha = alpha_init, Gamma2 = Gamma2_init,
  sigma2 = sigma2_init
)
# list of initialisations of the parameters to be estimated

hyperparam <- list(
  nu0 = nu0, nu1 = nu1, nu_Gamma = nu_Gamma, lb_Gamma = lb_Gamma, nu_sigma = nu_sigma,
  lb_sigma = lb_sigma, a = a, b = b, sigma2_mu = sigma2_mu, psi = psi, tau = tau
)
# list of hyperparameters

niter <- 500 # number of iterations of the algorithm
nburnin <- 350 # number of burn-in iterations (see Appendix A.1)
niterMH_phi <- 1 # number of iterations of the Metropolis Hastings algorithm for
# phi at S-step

source("scripts/Functions_SAEMVS.R")

res <- SAEM_MAP(niter, nburnin, niterMH_phi, Y, t, id, V_tilde, param_init, hyperparam,seed)

Delta <- 10^(seq(-2, 2, length.out = 20))

res <- Model_selection(Delta, niter, nburnin, niterMH_phi, Y, t, id, V_tilde, param_init, hyperparam,seed)

save(res,file=paste0("outputs/buildingResults/simulation/SAEMVS/res_",arr,".RData"))
