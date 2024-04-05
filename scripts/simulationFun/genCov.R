#' Generate a covTable following covSettings information.
#'
#' @param covSettings the covariates settings to add. The setting must be a list of covariates distribution information indexed by their names. The covariate distribution information is a list with two elements . The first one, named "distribution", must be a character string among c("beta","binom","cauchy","chisq","exp","f","gamma","geom","hyper","lnorm","nbinom","norm","pois","sample","unif","weibull"); the "sample" suffixes corresponds to the case where covariates are described by categorical covariates, eventualy not numerical, sampled according to a given set of probabilities.The rest of the available distributions correspond to the suffix of implemented distribution in the stats package. The second element of the list, named "elements", is the named parameters needed to call the sample function of the distribution (cf stats package). For the suffixe "sample", the needed arguments are "cat" corresponding to the covariate categories, and their respectives sample probabilities "prob". If given, the vector must be of same length, and if one (and only one) of them is set to NULL, they are replaced either with numerical categories from 1 to the length of prob or with equiprobabilities. (see examples). If a argument "bounds" is add the elements list, with either (or both) a min or a max value, the covariate will be sample in the $[min;max]$ set, with probability to be less or equal to min replaced by equal to min.
#' @param individualsNumber number of individuals to which simulate covariates.
#' @param covariateChar a single character to set to covariates unamed in the settings list. It will be concatenate with their number. (default to "X", usefull cause in addCov, the name of covariates are important to avoid repetition).
#'
#' @return a data.frame with "id" of individuals and associated covariates in columns.
#'
#' @examples
#' ## For stats function :
#'
#' ?runif #  runif(n, min = 0, max = 1) with "min and max : lower and upper limits of the distribution. Must be finite."
#' AGE = list(distribution = "unif" , elements =list(min=0,max=1))
#'
#' ?rnorm # rnorm(n, mean = 0, sd = 1) with "mean       : vector of means ; sd  : vector of standard deviations."
#' SIZE = list(distribution = "norm", elements = list(mean=165,sd=10), bounds=list(min=110,max=200)) #the size follow a normal distribution restrict between 110cm and 200cm.
#'
#' ## For sample :
#'
#' RACE = list(distribution = "sample", elements = list(cat=c("European","African"))) #prob set to c(0.5,0.5)
#' SMOKER = list(distribution = "sample", elements = list(cat=c("SMOKER","NON-SMOKER"),prob=c(0.3,0.7)))
#' STUDY = list(distribution ="sample", elements = list(prob=c(0.505,0.134,0.133,0.222,0.006)))  # cat set to c(1,2,3,4,5)
#'
#'
#' covSettings = list(AGE=AGE,SIZE=SIZE,RACE=RACE,SMOKER=SMOKER,STUDY=STUDY,WEIGHT = list(distribution="norm",elements=list(mean=70,sd=10),bounds=list(min=45)))
#'

############################################################################################################
genCov <- function(covSettings,individualsNumber,covariateChar="X"){
  #checkSettings(covSettings)
  covariateChar = paste0(covariateChar,1:length(covSettings))
  
  n=length(covSettings)
  covList = names(covSettings)
  
  covTable = data.frame(id=1:individualsNumber)
  
  for(i in 1:n){
    covInfo = covSettings[[i]]
    covName = covList[[i]] #set to NULL if nothing is provided and in this case name[i] will be taken
    if(is.null(covName)){covName <- covariateChar[i]}
    
    distribution = covInfo$distribution
    if(distribution=="uniform"){
      distribution="unif"
    }else if(distribution=="normal"){
        distribution="norm"
    }else if(distribution=="poisson"){
        distribution="pois"
      }
    elements = covInfo$elements
    bounds = covInfo$bounds
    
    
    
    if(distribution!="sample"){
      parameter = character(0)
      for(j in 1:length(elements)){
        if(j!=1) parameter = paste0(parameter,",")
        parameter = paste0(parameter,names(elements[j]),"=",elements[[j]])
      }
      
      command = paste0("r",distribution,"(",individualsNumber,",",parameter,")")
    }else{
      cat = elements$cat
      prob = elements$prob
      if(is.null(elements$cat) & is.null(elements$prob)){ stop(" Provide at leats one arguments among 'cat' and 'prob' for sample method.")
      }else if(is.null(cat)){ cat <- 1:length(prob)
      }else if(is.null(prob)){ prob <- rep(1/length(cat),length(cat))
      }else if(length(cat)!= length(prob)){ stop(" cat and prob must be of same length")}
      command = "sample(cat,individualsNumber,replace=TRUE,prob=prob)"
    }
    
    Cov = eval(parse(text=paste0(command)))
    
    if(!is.null(bounds)){
      lb = bounds$min
      ub = bounds$max
      if(!is.null(lb)) Cov[Cov<=lb] <- lb
      if(!is.null(ub)) Cov[Cov>=ub] <- ub
    }
    covTable <- cbind(covTable,Cov)
    colnames(covTable)[ncol(covTable)] <- covName
  }
  
  return(covTable)
}

checkSettings <- function(covSettings){
  if(typeof(covSettings)!="list"){stop("The covariates settings must be a list.")}
  flag=FALSE
  
  n = length(covSettings)
  listDist = unlist(lapply(covSettings,function(x){x$distribution}),use.names = FALSE)
  
  if(length(setdiff(listDist,c("beta","nbinom","hyper","binom","cauchy","chisq","exp","f","gamma","geom","norm","pois","sample","unif")))!=0){stop(paste0(" Unknown distribution provided : ",paste0(setdiff(listDist, c("beta","binom","cauchy","chisq","exp","f","gamma","geom","hyper","lnorm","nbinom","norm","pois","sample","unif")),collapse=", "),"."))}
  
  for(i in 1:n){
    cov = covSettings[[i]]
    distribution =cov$distribution
    elmts = names(cov$elements)
    if(is.null(elmts)){
      stop(" Please provided as indexed the names of the distribution parameters.")
    }else if(any(is.na(elmts))){stop(" Please provided as indexed the names of ALL the distribution parameters. Names of ",which(is.na(elmts)),"th parameters is missing.")}
    
    if(distribution=="sample"){
      rightElmts=c("cat","prob")
      if(length(elmts)==0){
        flag=TRUE
        warning(paste0(" Not enough parameters provided for ",i,"th covariates. For the ",distribution," distribution, parameters should be at least one among ", paste0(c("cat","prob"),collapse=" and "),"."))
      }else if(length(elmts)==1 && elmts!="cat" && elmts!="prob"){
        flag=TRUE
        warning(paste0("For ", distribution," distribution, parameters must be ",paste0(rightElmts,collapse=", ")," and not ",paste0(elmts,collapse=", "),"."))
      }else if(length(elmts)>2 && !(all(elmts %in% rightElmts) && all(rightElmts %in%elmts))){
        flag=TRUE
        warning(paste0("For ", distribution," distribution, parameters must be ",paste0(rightElmts,collapse=", ")," and not ",paste0(elmts,collapse=", "),"."))
      }
    }else{
      if(distribution == "beta"){
        rightElmts = c("shape1","shape2")
      }else if(distribution == "binom"){
        rightElmts = c("size","prob")
      }else if(distribution == "cauchy"){
        rightElmts =c("location","scale")
      }else if(distribution == "chisq"){
        rightElmts = c("df")
      }else if(distribution =="exp"){
        rightElmts = c("rate")
      }else if(distribution =="f"){
        rightElmts = c("df1","df2")
      }else if(distribution =="gamma"){
        rightElmts = c("shape","scale")
      }else if(distribution =="geom"){
        rightElmts = c("prob")
      }else if(distribution =="hyper"){
        rightElmts = c("m","n","k")
      }else if(distribution =="lnorm"){
        rightElmts = c("meanlog","sdlog")
      }else if(distribution =="nbinom"){
        rightElmts = c("size","prob")
      }else if(distribution =="norm"){
        rightElmts = c("mean","sd")
      }else if(distribution =="pois"){
        rightElmts = c("lambda")
      }else if(distribution =="unif"){
        rightElmts = c("min","max")
      }
      if(length(elmts)!=length(rightElmts)){stop(paste0(" Not enough parameters provided for ",i,"th covariates. For the ",distribution," distribution, parameters should be ", paste0(rightElmts,collapse=", "),"."))}
      
      if(!(all(elmts %in% rightElmts) && all(rightElmts %in%elmts))){
        flag=TRUE
        warning(paste0("For ", distribution," distribution, parameters must be ",paste0(rightElmts,collapse=", ")," and not ",paste0(elmts,collapse=", "),"."))
      }
    }
    if(flag){stop("Parameters elements for the covariate distribution aren't well-defined.")}
  }
}

genCorCov <- function(covSettings, corcovSettings, individualsNumber=NULL,oldCov=NULL, unamedCovariateChar="X"){
  if(!is.null(oldCov)){
    individualsNumber = nrow(oldCov)
    covTable = cbind(oldCov,genCov(covSettings,individualsNumber,unamedCovariateChar)[,-1])
  }else{
    if(is.null(individualsNumber)){stop("Please provide a number for individuals.")}
    covTable = genCov(covSettings,individualsNumber,unamedCovariateChar)
  }
  if(!is.null(corcovSettings)){
    covTable = data.table::setDT(covTable)
    
    n=length(corcovSettings)
    covList = names(corcovSettings)
    
    if(length(corcovSettings)>=1){
      covInfo = corcovSettings[[1]]
      distribution = covInfo$distribution
      formula = covInfo$formula
      variance = covInfo$variance
      link = covInfo$link
      
      if(is.null(variance)){variance = 0}
      if(is.null(link)){link = "identity"}
      if(is.null(distribution)){distribution = "normal"}
      
      def <- simstudy::defDataAdd(varname=covList[[1]],formula=formula,variance=variance,link=link,dist=distribution)
      
      if(length(corcovSettings)>1){
        for(i in 2:n){
          covInfo = corcovSettings[[i]]
          covName = covList[[i]]
          
          distribution = covInfo$distribution
          formula = covInfo$formula
          variance = covInfo$variance
          link = covInfo$link
          
          if(is.null(variance)){variance = 0}
          if(is.null(link)){link = "identity"}
          if(is.null(distribution)){distribution = "normal"}
          
          command = "def <- simstudy::defDataAdd(def,varname = covName,formula=formula,variance=variance,dist=distribution,link=link)"
          
          eval(parse(text=paste0(command))) 
        }
      }
      
      covTable = data.table::setDF(simstudy::addColumns(def, covTable))
    }
    
    boundedCov = names(which(lapply(lapply(corcovSettings,function(x){x$bounds}),is.null)==FALSE))
    for(c in boundedCov){
      min = corcovSettings[[c]]$bounds$min
      max = corcovSettings[[c]]$bounds$max
      if(!is.null(min)){
        covTable[covTable[[c]]<=min,c] <- min
      }
      if(!is.null(max)){
        covTable[covTable[[c]]>=max,c] <- min
      }
    }
  }
  return(covTable[,-which(colnames(covTable) %in% setdiff(colnames(oldCov),"id"))])
}