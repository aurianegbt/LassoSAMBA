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
