## Fun summarize and call
source("scripts/buildFS.R")

# Function just modified
source("scripts/modelBuildingFun/buildmlx.R")
source("scripts/modelBuildingFun/covariateModelSelection.R")
source("scripts/modelBuildingFun/covariateModelSelection.reg.R")

# Function implemented
source("scripts/modelBuildingFun/modelFromSelection.R")
source("scripts/modelBuildingFun/updateCov0.R")


source("scripts/modelBuildingFun/covariateModelSelection.lasso.R")
source("scripts/modelBuildingFun/applyMethodLasso.R")
source("scripts/modelBuildingFun/lassoSelectionPAR.R")

source("scripts/modelBuildingFun/covariateModelSelection.elasticnet.R")
source("scripts/modelBuildingFun/applyMethodElasticnet.R")
source("scripts/modelBuildingFun/elasticnetSelectionPAR.R")

p.weight <- Rsmlx:::p.weight
