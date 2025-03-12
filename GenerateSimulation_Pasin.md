# Generation of simulated datasets for Vaccinology framework with non gaussian correlated covariates (Pasin et al., 2019)

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
where $\eta^{\varphi_S}_i\sim^{iid}\mathcal N(0,\omega_{\varphi_S}^2)$, $\eta^L_i\sim^{iid}\mathcal N(0,\omega_L^2)$, $\eta^{Ab}_i\sim^{iid}\mathcal N(0,\omega_{Ab}^2)$. The observation are the defined as 
```math
Y_{ij} = \log_{10}(Ab_i(t_{ij}))+\varepsilon_{ij}
```
where $\varepsilon_i\sim^{iid}\mathcal N(0,\Sigma=\sigma^2_{Ab}I_{n_i})$.

The value used for the parameter are the one estimated from the EBOVAC trial (Eurosurveillance editorial team, 2014;Alexandre et al., 2023;Pasin et al., 2019). We then add up to 200 correlated covariates from a various range of distribution.


```r
set.seed(81511807)
dir <- function(path){if(!dir.exists(path)){dir.create(path)}}


## Load the library required 
suppressMessages({
  library(dplyr)
  library(lixoftConnectors)
  initializeLixoftConnectors("simulx")
  loadProject("data/simulationSetup/Pasin.smlx") # Simulx project with mechanistic model and parameters value 
  library(simstudy)
  library(data.table)
  source("scripts/simulationFun/genCov.R") 
  source("scripts/simulationFun/randomCovariate.R")
  load("data/simulationSetup/distribPasin.RData") # Covariance and correlation matrix based on real-data from Prevac-up clinical trial 
})

nb_ind <- 100
nb_replicates <- 10
nb_cov <- 200

```

For this scenario  with random distribution, we implemente a function \texttt{randomCovariate} providing a random distribution. If  This function outputs a list containing the name of the distribution and the necessary elements for sampling it in R~\citer. For example, if the distribution is ``unif'', for a uniform distribution, the elements will be ``min'' and ``max'' for the interval bounds. Table~\ref{tab:generation} summarizes the process of generating each of the distributions. To generate a distribution for a random variable $X$, the function call allows for a uniform draw from the possible distributions mentioned in table~\ref{tab:generation}, and then it draws the parameters according to the specified distributions, and the correlation matrix provided with the package \texttt{simstudy}~\cite{simstudy}.PK The covariates for each individual are consequently drawn according to this distribution.

| Distribution | Elements    | Generation                   |
|--------------|-------------|------------------------------|
| Gamma        | "shape"     | $\mathcal{P}(10)$            |
|              | "scale"     | $\mathcal{E}(0.5)$           |
| Normal       | "mean"      | $\mathcal{N}(0,1)$           |
|              | "sd"        | $\mathcal{E}(0.1)$           |
| Poisson      | "lambda"    | $\mathcal{E}(0.1)$           |
| Uniform      |             | $x = \mathcal{N}(0_2, I_2)$  |
|              | "min"       | $\min(x)$                    |
|              | "max"       | $\max(x)$                    |

**Table 1:** Process of generating various distributions using the implemented `randomCovariate` function.

Using this process, we generate the other covariate distribution. 

```r
## Generate distribution 
distribution = randomCovariate(n=nb_cov -3)


def <- defData(varname="AGE",formula =35,variance=16, dist = "normal")
def <- defData(def,varname="G1",formula=0,variance=1,dist="normal")
def <- defData(def,varname="G2",formula=0,variance=1,dist="normal")

for(i in 1:(nb_cov -3)){
  d= distribution[[i]]
  if(d$distribution=="gamma"){
    theta = d$elements$scale
    k = d$elements$shape
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=k*theta,variance = theta)
  }else if(d$distribution=="normal"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=d$elements$mean,variance=d$elements$sd)
  }else if(d$distribution=="poisson"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=d$elements$lambda)
  }else if(d$distribution=="uniform"){
    def <- defData(def,varname=paste0("Gen",i),dist=d$distribution,
                   formula=paste0(d$elements,collapse = ";"))
  }
}
setNbReplicates(1)
genCorMat <- genCorMat + 1e-10 * diag(ncol(genCorMat))
covTableALL = genCorFlex(nb_replicates*nb_ind,def,corMatrix=genCorMat)
``` 

The simulated covariates are then add to the simulX project to generate the individual dynamics. 

```r
for(i in 1:nb_replicates){
  covTable = covTableALL[(1+(i-1)*nb_ind):(i*nb_ind),] %>% mutate(id = (id-1)%%nb_ind+1)

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
  
  dir("data/simulationFiles/FilesPasin")
  dir("data/simulationFiles/FilesPasin/covTable") # used for graphs 
  dir("data/simulationFiles/FilesPasin/simulation")
  
  if(i ==1){
    write.csv(covTable,file=paste0("data/simulationFiles/FilesPasin/covTable/covTable_",i,".txt"),quote = F,row.names = F)
  }
  write.csv(dataset,file=paste0("data/simulationFiles/FilesPasin/simulation/simulation_",i,".txt"),quote = F,row.names = F)
  
  unlink(paste0("tmpfile",i,".txt"))
}

headerTypes = c("id","time","observation",rep("contcov",3),sapply(sapply(distribution,FUN=function(x){x$distribution})=="poisson",FUN=function(x){if(x){"catcov"}else{"contcov"}}))
save(headerTypes,file="data/simulationFiles/FilesPasin/headerTypes.RData")
```
