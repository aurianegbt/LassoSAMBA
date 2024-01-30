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
source("Fun/Rsmlx/covariateModelSelection.rlasso.R")
source("Fun/Rsmlx/applyMethodLasso.R")
source("Fun/Rsmlx/lassoSelectionPAR.R")

source("Fun/Rsmlx/covariateModelSelection.elasticnet.R")
source("Fun/Rsmlx/covariateModelSelection.relasticnet.R")
source("Fun/Rsmlx/applyMethodElasticnet.R")
source("Fun/Rsmlx/elasticnetSelectionPAR.R")

source("Fun/Rsmlx/covariateModelSelection.sharp.R")
source("Fun/Rsmlx/applyMethodsharp.R")

source("Fun/Rsmlx/covariateModelSelection.rsharp.R")
source("Fun/Rsmlx/applyMethodrsharp.R")
source("Fun/Rsmlx/sharpModified/SerialRegression.R")
source("Fun/Rsmlx/sharpModified/VariableSelection.R")


p.weight <- Rsmlx:::p.weight
PenalisedRegression <- sharp::PenalisedRegression
