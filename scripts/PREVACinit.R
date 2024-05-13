library(lixoftConnectors)
initializeLixoftConnectors(force=TRUE)
library(dplyr)

load(file="data/applicationFiles/arm1/antibodyData.RData")

removeID = unique(database[is.na(database$init),"id"])

database = database %>% 
  mutate(time=stringr::str_remove_all(time," days")) %>%
  filter(!(id %in% removeID))%>%
  filter(!is.na(yAb)) %>%
  filter(id != 30402119)

write.csv(database,"data/applicationFiles/antibodyArm1.txt",quote = F,row.names = F)

newProject(modelFile="data/modelFiles/PasinReg.txt",
           data=list(dataFile="data/applicationFiles/antibodyArm1.txt",
                     headerTypes = c("id","ignore","time","observation","regressor")))

setIndividualParameterVariability(delta_L = FALSE)
setPopulationParameterInformation(delta_L_pop = list(initialValue=0.000316,method="FIXED"),
                                  delta_S_pop = list(initialValue=0.231,method="FIXED"))

setErrorModel(yAb="constant")

scenario <- getScenario()
scenario$tasks["standardErrorEstimation"] <- TRUE 
setScenario(scenario)

setPopulationParameterInformation(getFixedEffectsByAutoInit())

runScenario()

setInitialEstimatesToLastEstimates()
PopulationParameterInformation <- getPopulationParameterInformation()
save(PopulationParameterInformation,file="data/applicationFiles/PopulationParameterInit.RData")
