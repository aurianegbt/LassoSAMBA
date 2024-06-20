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

data.dir<- "data/simulationFiles/FilesPKcorr/simulation"
dir(data.dir)
data.dir.cov <- "data/simulationFiles/FilesPKcorr/covTable"
dir(data.dir.cov)
project_ini.dir <- "project_ini"
project_star.dir <- "project_star"
file.results <- "simulations.RData"

#Simulation settings
M <- 100 #number replicates
seed0 <- 12345

N <- 100 #Number of individuals
ncov <- 200 #Number of covariates
# Correlation matrix
# load("data/applicationFiles/arm1/DataTransCoding.RData")
# 
# aux = (dataTransCoding %>% filter(visit=="D63"))[,c(6,9:16894)]
# genesKept = sample(names(sort(apply(aux[,-1],2,sd),decreasing = TRUE))[2:201]) # keep only 200
# 
# aux <- aux[,genesKept]
# colnames(aux) <- c(1:ncol(aux))
# 
# genCorMat <- cor(aux,method="spearman")
# mu = c(0,0,rowMeans(t(aux[,-c(1:2)]))) # C1, C2 have already means defined in simulX
# mu = rep(0,ncol(genCorMat))
# genCovMat = sd %*% genCorMat %*% sd
# 
# epsilon <- 1e-10
# genCorMat <- genCorMat + epsilon * diag(ncol(genCorMat))
# 
# save(mu,genCorMat,genCovMat,file="data/simulationSetup/distribPKcorr.RData")
# 
# 
# ## Plot
# corrplot =  ggcorrplot(genCorMat, colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used.",
#                                 face="bold",size=50,color="#882255"))
# ggsave("outputs/figures/explanatory/corrPKcorr.png",
#        height = 10000, width = 10000, units = "px", bg='transparent')
# 
# corrplot =  ggcorrplot(genCorMat[1:50,1:50], colors= c("#446494","#eeeeee","#882255"))  + theme(plot.title = element_text(size=20, face="plain")) + theme(plot.title = element_text(size=20, face="plain"))
# 
# annotate_figure(corrplot,
#                 top = text_grob("Theoretical Correlation Matrix used, zommed on the 50 first covariates.",
#                                 face="bold",size=50,color="#882255"))
# ggsave("outputs/figures/explanatory/corrPKcorrzoom.png",
#        height = 10000, width = 10000, units = "px", bg='transparent')

load(file="data/simulationSetup/distribPKcorr.RData")

## Generate 200 correlated covariates and then create 100 replicates
t.y <- c(0.25, 0.5, 1, 2, 5, 8, 12, 16, 20, 24, 30) #observation times
D <- 1000 

# Model specification
lib.1cpt <- 'lib:oral1_1cpt_kaVCl.txt'
p.1cpt <- c(ka_pop=1, V_pop=10, Cl_pop=2,
            omega_ka=0.2, omega_V=0.3, omega_Cl=0.3,
            corr_V_Cl=0.6,
            beta_V_C1=0.2, beta_Cl_C1=-0.2, beta_Cl_C2=0.3,
            a=2, b=0.1)

PKmodel.1cpt <- inlineModel("
[COVARIATE]
input = {C1, C2}

[INDIVIDUAL]
input = {ka_pop, omega_ka, Cl_pop, omega_Cl, V_pop, omega_V, 
C1, C2, beta_V_C1, beta_Cl_C1, beta_Cl_C2, corr_V_Cl}

DEFINITION:
ka = {distribution=logNormal, typical=ka_pop, sd=omega_ka}
V  = {distribution=logNormal, typical=V_pop, covariate=C1, coefficient=beta_V_C1, sd=omega_V}
Cl = {distribution=logNormal, typical=Cl_pop, covariate={C1, C2}, coefficient={beta_Cl_C1, beta_Cl_C2}, sd=omega_Cl}
correlation = {level=id, r(V, Cl)=corr_V_Cl}

[LONGITUDINAL]
input = {ka, V, Cl, a, b}

EQUATION:
Cc = pkmodel(ka, V, Cl)
g = sqrt(a^2 + (b*Cc)^2)

DEFINITION:
y = {distribution=normal, prediction=Cc, sd=g}
")

p <- p.1cpt
p.pop <- p[1:3]
p.name <- c("ka", "V", "Cl")
structural.model <- lib.1cpt
PKmodel <- PKmodel.1cpt

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
  covariate <- data.frame(id=factor(1:N), mvtnorm::rmvnorm(N,mean=mu,sigma=genCorMat))
  names(covariate)[-1] <- cov.name
  if(m==1){
    write.csv(covariate, file=paste0(data.dir.cov,"/covTable_1.txt"),quote = F,row.names = F)
  }
  
  
  r <- simulx(model     = PKmodel, 
              parameter = p,
              covariate = covariate[c('id','C1', 'C2')],
              treatment = list(time=0,  amount=D),
              output    = list(name="y", time=t.y),
              settings  = list(seed=seed[m]))
  
  d <- bind_rows(r$y,r$treatment %>% mutate(id=as.factor(id)) %>% select(c(id, time,amount)))
  d <- left_join(d, covariate, by="id")
  d[is.na(d)] <- "."
  write.csv(d, file=file.path(data.dir,data.files[m]), quote=F, row.names=F)
}
