library(lixoftConnectors)
library(sharp)
library(foreach)
library(dplyr)
source("scripts/modelBuildingFun/sharpCalibrationPlot.R")
initializeLixoftConnectors()

sim = read.csv(file="data/simulationFiles/FilesPKcorr/simulation/simulation_1.txt") %>% select(id,time,y,amount)%>% filter(id<=95)

load("data/applicationFiles/DataTransCodingD63.RData")

dir.create("/beegfs/agabaut/tmpNOISYSIM/")
dataSet = merge(sim,genesCodingD63[,1:300] %>% mutate(id=1:94),by = "id")

write.csv(dataSet,file = "/beegfs/agabaut/tmpNOISYSIM/dataset.txt",quote = F,row.names = F)

newProject(modelFile = "data/modelFiles/PKcorr.txt",
           data=list(dataFile="/beegfs/agabaut/tmpNOISYSIM/dataset.txt",
                     headerTypes=c("id","time","observation","amount",rep("contcov",299))))

saveProject(projectFile = "/beegfs/agabaut/tmpNOISYSIM/BUILD.mlxtran")

dir.create("outputs/figures/explanatory/CalibrationPlot/Noisy/")
project = "/beegfs/agabaut/tmpNOISYSIM/BUILD.mlxtran"
pathToSave = "outputs/figures/explanatory/CalibrationPlot/Noisy/"
title.param =c(Cl = latex2exp::TeX(r"(Calibration plot for $Cl$ parameters)"),
               ka= latex2exp::TeX(r"(Calibration plot for $ka$ parameters)"),
               V = latex2exp::TeX(r"(Calibration plot for $V$ parameters)"))

pathToSave = "outputs/figures/explanatory/CalibrationPlot/Noisy/FDP5"
sharpCalibrationPlot.init(project,pathToSave,title.param,JPEG = T,PNG=F,FDP_thr=0.05)

pathToSave = "outputs/figures/explanatory/CalibrationPlot/Noisy/FDP50"
sharpCalibrationPlot.init(project,pathToSave,title.param,JPEG = T,PNG=F,FDP_thr=0.50)

pathToSave = "outputs/figures/explanatory/CalibrationPlot/Noisy/FDP80"
sharpCalibrationPlot.init(project,pathToSave,title.param,JPEG = T,PNG=F,FDP_thr=0.80)

