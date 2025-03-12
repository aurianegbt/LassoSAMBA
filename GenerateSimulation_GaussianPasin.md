# Generation of simulated datasets for Vaccinology framework with gaussian correlated covariates (Pasin et al., 2019)

We simulate for 100 individuals the antibody production by considering two Antibodies secreting cells (ASC), denoted by S -\textit{for short-live}- and L -\textit{for long-live}- (at rates $\varphi_S$ and $\varphi_L$ resp.) and characterized by their half-life ($\delta_S$ and $\delta_L$ resp.). Antibodies are supposed to decay at rate $\delta_{Ab}$. We add significant covariates on $\varphi_S$, $\varphi_L$ and $\delta_{Ab}$ parameters, and then add noisy genes in order to have finally 200 covariates. We create two scnearios, one with only standard gaussian covariates and one with randomly choose distribution, with a mix of categorical, and continuous covariates for sake of validity. For these two scenarios, the covariates are correlated. The generation process for these covariates is detailed later. The mechanistic model is then : 
```math
\forall i\leq N,j\leq n_i,   \left\{\begin{array}{rcl}
    \frac{d}{dt} Ab_i(t_{ij}) &=& {\varphi_S}_i e^{-\delta_S t_{ij}} + {\varphi_L}_i e^{-\delta_L t_{ij}} - {\delta_{Ab}}_i Ab_i(t_{ij}) \\
    Ab_i(t_{i0}=0) &=& {Ab_0}
\end{array}\right.
```
with 
```math
\displaystyle\left\{
\begin{array}{rcl}
         \log({\varphi_S}_i) &=& \log({\varphi_S}_{pop}) + \eta^{\varphi_S}_i \\
         \log({\varphi_L}_i) &=& \log({\varphi_L}_{pop})  + \eta^L_i \\
         \log({\delta_{Ab}}_i) &=& \log({\delta_{Ab}}_{pop})   +\eta^{Ab}_i
    \end{array}\right.
```
where $\eta_i^{\varphi_S}\sim^{iid}\mathcal{N}(0,\omega^2_{\varphi_S})$, $\eta^L_i\sim^{iid}\mathcal{N}(0,\omega_L^2)$, $\eta_i^{Ab}\sim^{iid}\mathcal{N}(0,\omega^2_{Ab})$. The observation are the defined as 
```math
Y_{ij} = \log_{10}(Ab_i(t_{ij}))+\varepsilon_{ij}
```
where $\varepsilon_i\sim^{iid}\mathcal N(0,\Sigma=\sigma^2_{Ab}I_{n_i}) $.

The value used for the parameter are the one estimated from the EBOVAC trial (Eurosurveil-lance editorial team, 2014;Alexandre et al., 2023;Pasin et al., 2019). We then add up to 200 correlated covariates.


```r
set.seed(1710)
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

nb_replicates <- 10
nb_ind <- 100
nb_cov <- 200
```

The correlation matrix used for the vaccinology simulation is based on correlation estimation of the genetic expression of Prevac-UP immuno substudy participants (Badio et al., 2021). We estimated correlation from the clinical trial data using the Spearman method (Spearman, 1904). This allows to have realistic correlation among our simulated genes. The following file contains covariance, correlation matrix and the mean value vector of covariates. We can also load the Simulx project created, with parameters values and mechanistic model already set. 

```r
load("data/simulationSetup/distribPasin.RData")
loadProject("data/simulationSetup/Pasin.smlx")
```

The covariates are simulated for the number of replicates desired and add to the simulx project.

```r
covTableALL = as.data.frame(mvtnorm::rmvnorm(n=nb_replicates*nb_ind,mean=mu,sigma = genCovMat))
colnames(covTableALL) <- c("AGE","G1","G2",paste0("Gen",1:(nb_cov-3)))

for(i in 1:nb_replicates){
  covTable = covTableALL[(1+(i-1)*nb_ind):(i*nb_ind),]
  covTable <- cbind(id=1:nb_ind,covTable)
  
  write.csv(covTable[,1:4],paste0("tmpfile",i,".txt"),quote = F,row.names = F)
  
  defineCovariateElement(name=paste0("covTable",i),
                         element = paste0("tmpfile",i,".txt"))
  
  if(i==1){
    setGroupElement(group=paste0("simulationGroup",i), elements = c(paste0("covTable",i)))
  }else{
    addGroup(paste0("simulationGroup",i))
    setGroupElement(group=paste0("simulationGroup",i), elements = c(paste0("covTable",i)))
  }
}
  
runSimulation()
sim <- getSimulationResults()
```

Once the simulation launch, data files are saved. As it is usually done before analyzis, the age of participants is centered.

```r
for(i in 1:nb_replicates){
  dataset = sim$res$yAB[sim$res$yAB$group==paste0("simulationGroup",i),c("original_id","time","yAB")] %>%
    rename(id=original_id)
  
  covTable = covTableALL[(1+(i-1)*nb_ind):(i*nb_ind),] %>% 
    mutate(AGE = AGE - mean(covTableALL[(1+(i-1)*nb_ind):(i*nb_ind),]$AGE)) %>%
    mutate(id=1:nb_ind,.before = AGE) %>%
    rename(cAGE = AGE)
  
  dataset = merge(dataset,covTable,by = "id")
  
  dataset$id <- as.numeric(dataset$id)
  
  dataset = dataset %>% arrange(id,time)
  
  dir("data/simulationFiles/FilesGaussianPasin")
  dir("data/simulationFiles/FilesGaussianPasin/covTable") # used for graphs so only need to save for first replicates 
  dir("data/simulationFiles/FilesGaussianPasin/simulation")
  
  if(i==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesGaussianPasin/covTable/covTable_",i,".txt"),quote = F,row.names = F) 
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesGaussianPasin/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}


headerTypes = c("id","time","observation",rep("contcov",nb_cov))
save(headerTypes,file="data/simulationFiles/FilesGaussianPasin/headerTypes.RData")
```
