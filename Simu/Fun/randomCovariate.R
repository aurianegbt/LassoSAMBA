randomCovariate <-function(n=1,NORM=FALSE){
  if( n %% 1 != 0 | n==0 ){stop("The number of covariates to generate must be an integer.")}
  res = list()
  for(i in 1:n){ res <- append(res,list(covRD_one(NORM)))}
  return(res)
}

covRD_one <- function(NORM=FALSE){
  if(NORM){
    mean = rnorm(1,sd=5)
    sd=rexp(1,0.2)
    elements=list(mean=mean,sd=sd)
    return(list(distribution="norm",elements=elements))
  }else{
    listDist = c("beta","binom","cauchy","chisq","exp","f","gamma","geom","norm","pois","sample","unif")
    distribution = sample(listDist,size=1)
    if(distribution == "beta"){
      shape1 = rpois(1,20)
      shape2 = rpois(1,20)
      elements=list(shape1=shape1,shape2=shape2)
    }else if(distribution == "binom"){
      prob = runif(1,min=0.2,max=0.8)
      size = max(min(rpois(1,10),20),1)
      elements=list(size=size,prob=prob)
    }else if(distribution == "cauchy"){
      location = rnorm(1,sd = 20)
      scale = rexp(1,0.1)
      elements=list(location=location,scale=scale)
    }else if(distribution == "chisq"){
      df = max(rpois(1,20),3)
      elements=list(df=df)
    }else if(distribution =="exp"){
      rate = rexp(1,2)
      elements=list(rate=rate)
    }else if(distribution =="f"){
      df1 = max(rpois(1,20),3)
      df2 = max(rpois(1,20),3)
      elements=list(df1=df1,df2=df2)
    }else if(distribution =="gamma"){
      shape = rpois(1,10)
      scale = rexp(1,1/2)
      elements=list(shape=shape,scale=scale)
    }else if(distribution =="geom"){
      prob = runif(1,max=0.6)
      elements=list(prob=prob)
    }else if(distribution =="sample"){
      size = max(rpois(1,5),2)
      prob = rep(0,size)
      for( i in 1:(size-1)){
        prob[i] <- runif(1,min=0,max=1-sum(prob))
      }
      prob[size]=1-sum(prob)
      elements=list(prob=prob)
    }else if(distribution =="norm"){
      mean = rnorm(1)
      sd=rexp(1,0.3)
      elements=list(mean=mean,sd=sd)
    }else if(distribution =="pois"){
      lambda=rexp(1,0.1)
      elements=list(lambda=lambda)
    }else if(distribution =="unif"){
      bounds = rnorm(2,sd = 10)
      min = round(min(bounds),digits = 2)
      max= round(max(bounds),digits = 2)
      while(max-min<0.3){
        bounds = rnorm(2,sd = 10)
        min = round(min(bounds),digits = 2)
        max= round(max(bounds),digits = 2)
      }
      elements=list(min=min,max=max)
    }
    return(list(distribution=distribution,elements=elements))
  }
}
