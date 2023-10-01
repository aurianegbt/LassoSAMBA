createDatasets <- function(pathToCov,
                           modelFile, 
                           outputElement, 
                           populationElement){
  
  
  coTable = read.csv(pathToCov)
  individualsNumber = as.integer(nrow(covTable))
  
  
  newProject(modelFile = modelFile)
  
  definePopulationElement(name=populationElement$name, 
                          element = populationElement$element)
  
  defineOutputElement(name=outputElement$name, 
                      element = outputElement$element)
  removeGroup("simulationGroup1")
  
  defineCovariateElement(name=paste0("covTable",i),
                         element =  pathToCov)
  
  grp=paste0("Simulation",sep="_",i)
  addGroup(grp)
  setGroupSize(grp,individualsNumber)
  setGroupElement(grp,elements = populationElement$name)
  setGroupElement(grp,elements=paste0("covTable",i))
  setGroupElement(grp,elements=outputElement$name)
  
  runSimulation()
  sim=getSimulationResults()
  y = outputElement$element$output
  res = sim$res[[y]]
  data_sim = merge(res,covTable,by="id")
  return(data_sim)
}