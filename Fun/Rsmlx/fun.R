## Fun summarize and call
source("Fun/buildFS.R")

# Function just modified
source("Fun/Rsmlx/buildmlx.R")
source("Fun/Rsmlx/covariateModelSelection.R")
source("Fun/Rsmlx/covariateModelSelection.reg.R")

# Function implemented

source("Fun/Rsmlx/modelFromSelection.R")
source("Fun/Rsmlx/updateCov0.R")


source("Fun/Rsmlx/covariateModelSelection.lasso.R")
source("Fun/Rsmlx/applyMethodLasso.R")
source("Fun/Rsmlx/lassoSelectionPAR.R")

source("Fun/Rsmlx/covariateModelSelection.StepAIC.R")
source("Fun/Rsmlx/applyMethodStepAIC.R")
source("Fun/Rsmlx/StepAICSelection.R")


source("Fun/Rsmlx/covariateModelSelection.elasticnet.R")
source("Fun/Rsmlx/applyMethodElasticNet.R")
source("Fun/Rsmlx/elasticnetSelectionPAR.R")

p.weight <- Rsmlx:::p.weight
