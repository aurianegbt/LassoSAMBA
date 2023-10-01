addCov <- function(pathToCov=NULL,pathToSaveCov=NULL,covNew,pathToSimS=NULL,pathToSaveSimS=NULL,print=TRUE){
  if(print){cat("- - - New covariate(s) initialization  - - -\n")}
  
  ## Handling null arguments : 
  if(is.null(pathToSaveSimS)){pathToSaveSimS = pathToSimS}
  if(!is.null(pathToCov)){
    covTable = read.csv(pathToCov)
    Nbr_ind = length(unique(covTable$id))
    if(nrow(covNew)<Nbr_ind){stop("Please provided enough new covariates to add. The number of individuals provided is greater than the covariate informations. ")}
    
    covToAdd = setdiff(colnames(covNew),colnames(covTable))
    
    if(length(covToAdd)==0){stop("No covariates provided to add.")}
    if(length(c("id",covToAdd))!=length(colnames(covNew))){warning("Some Covariates in covNew are already defined in covTable. They will be ignored.")}
    
    if(is.null(pathToSaveCov)){pathToSaveCov=pathToCov}
    covTable_new = covNew[covNew$id <= Nbr_ind,c("id",covToAdd)]
    
    
    covTable <- merge(covTable,covTable_new,by="id")
  }else{
    covTable_new = covNew
    covTable = covTable_new
    covToAdd = setdiff(colnames(covTable),"id")
  }
  if(!is.null(pathToSaveCov)){write.csv(covTable,file=pathToSaveCov,quote=FALSE,row.names=FALSE)}
  
  
  if(print){
    if(length(covToAdd)==1){cat(paste0("One covariate to add :\n",covToAdd,".\n"))
   }else{
        cat(paste0(length(covToAdd)," covariates to add :\n",paste0(covToAdd[-length(covToAdd)],collapse=", ")," and ",covToAdd[length(covToAdd)],".\n"))}
  }
  if(print){cat("- - - Starting covariate concatenation - - -\n")}
  
  ## Function core : 
  if(!is.null(pathToSimS)){
    for (i in 1:length(pathToSimS)){
      path = pathToSimS[i]
      pathToSave = pathToSaveSimS[i]
      if(print){cat(paste0("Work on ",path, "...\n"))}
      data_sim = read.csv(path)
      
      data_sim <- merge(data_sim,covTable_new,by="id")

      write.csv(data_sim,file=pathToSave,quote=FALSE,row.names=FALSE)
    }
  }
  cat(" - - - - - - - - - - - - - - - - - - - - - -\n")
}