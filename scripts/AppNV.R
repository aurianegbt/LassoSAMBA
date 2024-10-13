# Empty memory ------------------------------------------------------------
rm(list=ls())
dir <- function(x){if(!dir.exists(x)){dir.create(x)}}

# Load packages, arguments, functions and  data  --------------------------
suppressWarnings({
  suppressMessages({
    library(ggplot2)
    library(cowplot)
    library(stringr)
    library(dplyr)
    library(mvnfast)
    library(doParallel)
    library(lixoftConnectors)
    library(foreach)
    library(dplyr)
    initializeLixoftConnectors(software = "monolix",
                               path = "/cm/shared/dev/modules/generic/apps/tools/monolix/2023R1bis/",
                               force = TRUE)
  })
})

arr <- as.numeric(slurmR::Slurm_env(x='SLURM_ARRAY_TASK_ID'))

args = commandArgs(trailingOnly=TRUE)
# args=c("reg","NULL","FALSE","'covariate'")
buildMethod = args[1]
FDP_thr = if(args[2]!="NULL"){as.numeric(args[2])}else{NULL}
noCov0 = as.logical(args[3])
eval(parse(text=paste0("model= ",args[4])))
ch = paste0("chr",paste0(1:7,rep(c("A","B","D"),each=7)))[arr]

source("scripts/MBFun.R")

# Presentation of Job -----------------------------------------------------
cat("===========================================================================\n
- - - - - - - - - - - - - - Work Information - - - - - - - - - - - - - - -\n")
cat("Application to Plant data\n")
cat(" • Chromosome : ",ch,"\n")
cat(" • Method used : ",buildMethod,"\n")
if(buildMethod=="sharp"){cat("     > FDP_thr=",FDP_thr*100,"% \n")}
cat(" • Model(s) built : ",paste0(model,collapse=", "))
if(noCov0){cat("\n Model building done without statistical tests.")}
cat("\n===========================================================================\n")

# Folder and paths  -------------------------------------------------------
pathToResults = paste0("outputs/buildingResults/application/ResultsPlant_",ch)
dir(pathToResults)
pathToResults = paste0(pathToResults,"/Results_",buildMethod,if(noCov0){"noCov0"},if(buildMethod=="sharp"){FDP_thr*100})
dir(pathToResults)
pathToResults = paste0(pathToResults,"/",paste0(model,collapse="_"),".RData")

# Temporary Files ---------------------------------------------------------
temporaryDirectory = paste0("/beegfs/agabaut/tmpPLANT_",ch,"_",buildMethod,if(noCov0){"noCov0"},if(buildMethod=="sharp"){FDP_thr*100})
if(dir.exists(temporaryDirectory)){
  warning("A folder already exists for this job. It has been deleted.")
  unlink(temporaryDirectory,recursive=T,force=T)
}

dir.create(temporaryDirectory)

if(!file.exists(paste0("data/applicationFiles/plant/chromosome/data_file",ch,".txt"))){
  
  load("data/applicationFiles/plant/data_obs.Rdata")
  load("data/applicationFiles/plant/data_ch.RData")
  #### Firt we construct the matrix V_tilde associated with the chromosome 6A:
  
  markers_ch=data_chromosome$name[(data_chromosome$V1==ch)]
  V_ch=SNPs[,colnames(SNPs)%in%markers_ch]
  n=dim(V_ch)[1]
  markers_noQTL=setdiff(colnames(SNPs),heading_QTL) #list of SNPs names that are not heading QTLs
  
  V_ch_QTL=as.matrix(V_ch[,colnames(V_ch)%in%heading_QTL])
  nb_QTL=dim(V_ch_QTL)[2] #number of heading QTLs in chr6A
  if (nb_QTL!=1){
    couples1 <- c()
    couples2 <- c()
    similarities <- matrix(NA,nb_QTL,nb_QTL)
    similarities2 <- matrix(NA,nb_QTL,nb_QTL)
    for (i in 1:(nb_QTL-1)){
      for (j in (i+1):nb_QTL){
        similarities[i,j] <- sum(V_ch_QTL[,i]==(V_ch_QTL[,j]))
        similarities2[i,j] <- sum(V_ch_QTL[,i]==(1-V_ch_QTL[,j]))
        if (similarities[i,j] == n){
          couples1 <- c(couples1,i,j)
        }
        if (similarities2[i,j] == n){
          couples2 <- c(couples2,i,j)
        }
      }
    }
    indremove1=c()
    if (is.null(couples1)==FALSE){
      couples1 <- matrix(couples1,ncol=2,byrow=T)
      v1 <- unique(couples1[,1])
      v2 <- unique(couples1[,2])
      indkeep1 <- setdiff(v1,v2)
      indremove1 <- v2
    }
    
    indremove2=c()
    if (is.null(couples2)==FALSE){
      couples2 <- matrix(couples2,ncol=2,byrow=T)
      w1 <- unique(couples2[,1])
      w2 <- unique(couples2[,2])
      indkeep2 <- setdiff(w1,w2)
      indremove2 <- w2
    }
    
    remove <- union(indremove1,indremove2)
    if (is.null(remove)==FALSE){
      keep=colnames(V_ch_QTL)[-remove]
      V_ch=V_ch[,colnames(V_ch)%in%c(markers_noQTL,keep)]
      markers_ch=colnames(V_ch)
    }
    
  }
  
  
  nb_QTL=sum(markers_ch%in%heading_QTL)
  
  ## And then for the others:
  couples1 <- c()
  nb <- dim(V_ch)[2]
  similarities <- matrix(NA,nb,nb)
  for (i in 1:(nb-1)){
    for (j in (i+1):nb){
      similarities[i,j] <- sum(V_ch[,i]==(V_ch[,j]))
      if (similarities[i,j] == n){
        couples1 <- c(couples1,i,j)
      }
    }
  }
  
  couples1 <- matrix(couples1,ncol=2,byrow=T)
  
  v1 <- unique(couples1[,1])
  v2 <- unique(couples1[,2])
  indkeep1 <- setdiff(v1,v2)
  #We check that we are not deleting any heading QTLs:
  keep=intersect(which(markers_ch%in%heading_QTL),v2)
  # keep=integer(0): it's ok, we are not deleting heading QTL
  indremove1 <- v2
  
  couples2 <- c()
  nb <- dim(V_ch)[2]
  similarities <- matrix(NA,nb,nb)
  for (i in 1:(nb-1)){
    for (j in (i+1):nb){
      similarities[i,j] <- sum(V_ch[,i]==(1-V_ch[,j]))
      if (similarities[i,j] == n){
        couples2 <- c(couples2,i,j)
      }
    }
  }
  
  couples2 <- matrix(couples2,ncol=2,byrow=T)
  w1 <- unique(couples2[,1])
  w2 <- unique(couples2[,2])
  #We check that we are not deleting any heading QTLs:
  keep=intersect(which(markers_ch%in%heading_QTL),w2)
  # keep=integer(0): it's ok, we are not deleting heading QTL
  indremove2 <- w2
  
  remove <- union(indremove1,indremove2)
  
  V_ch <- V_ch[,-remove]
  nb_QTL=sum(heading_QTL%in%colnames(V_ch))
  ind_QTL=which(colnames(V_ch)%in%heading_QTL) #indices of the columns of V_ch that correspond to heading QTLs
  dim(V_ch)
  
  sspop <- readRDS("data/applicationFiles/plant/resPCOAdf.Rds") #5 new covariates to control for subpopulation structure, obtained by PCA
  V_tilde_ch= as.matrix(cbind(Intercept=rep(1,n),sspop,V_ch))
  
  V_tilde_ch[,-1] <- scale(V_tilde_ch[,-1],center=T,scale=T)
  data_ch=list(V_tilde=V_tilde_ch,nb_QTL=nb_QTL,ind_QTL=ind_QTL+6,markers=colnames(V_ch))
  
  ## creation of the dataset 
  
  varieties=unique(data_obs$GENOTYPE)
  obs_date=unique(data_obs$Day)
  n=length(unique(data_obs$GENOTYPE))
  J=length(obs_date)
  Id=rep(c(1:n),each=J)
  id=as.matrix(Id)
  Y=c()
  for (i in varieties){
    for (j in obs_date){
      Y=c(Y,data_obs[(data_obs$GENOTYPE==i)&(data_obs$Day==j),]$Y)
    }
  }
  t=rep(unique(data_obs$Day),n)
  
  V_tilde=cbind(id=1:n,data_ch$V_tilde[,-1])
  
  data = merge(cbind(id=Id,time=t,obs=Y),V_tilde,by="id")
  write.csv(data,file=paste0("data/applicationFiles/plant/chromosome/data_file_",ch,".txt"),quote = F,row.names = F)
}

sim = read.csv(paste0("data/applicationFiles/plant/chromosome/data_file_",ch,".txt"))
nbCov = ncol(sim) - 3

newProject(modelFile = "data/modelFiles/NAVEAU.txt",
           data = list(dataFile=paste0("data/applicationFiles/plant/chromosome/data_file",ch,".txt"),
                       headerTypes = c("id","time","observation",rep("contcov",nbCov))))


setIndividualParameterVariability(psi1=FALSE)
setPopulationParameterInformation(psi1_pop=c(initialValue=100,method="FIXED"))
setIndividualParameterDistribution(phi="normal")
setIndividualParameterDistribution(psi2="normal")
setErrorModel(yAb="constant")

paramToUse = "phi"
saveProject(paste0(temporaryDirectory,"/Build.mlxtran"))

# Model Building
res = buildmlx(project = paste0(temporaryDirectory,"/Build.mlxtran"),
               buildMethod = buildMethod,
               model=model,
               test=FALSE,
               thresholdsSS=thresholdsSS,
               FDP_thr = FDP_thr,
               p.max=if(noCov0){1}else{0.1})

Model <- Rsmlx:::mlx.getIndividualParameterModel()
covModel <- lapply(Model$covariateModel,FUN=function(cov){names(cov[which(cov)])})
time=res$time
iter=res$iter
LL = lixoftConnectors::getEstimatedLogLikelihood()

save(res,Model,covModel,time,iter,LL,file=pathToResults)
