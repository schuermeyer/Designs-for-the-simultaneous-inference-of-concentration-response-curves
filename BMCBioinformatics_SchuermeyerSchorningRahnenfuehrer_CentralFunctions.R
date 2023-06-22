# First the central functions are introduced: 


# Fisherinformation matrix depending on w,x, and theta
info <- function(w,x,gradient,...){
  dimTheta<-length(gradient(1,...))
  tmp<-matrix(0,dimTheta,dimTheta)
  for(i in 1:length(x)){
    tmp=tmp+w[i]*(gradient(x[i],...)%*%t(gradient(x[i],...)))
  }
  return(tmp)
}

# D-optimality criterion
Dcrit<-function(w,x,gradient,...){
  return((det(info(w,x,gradient,...)))^(1/4))
}

# Gradient of sigmoid Emax from DoseFinding package
mygrad = function(x, theta)
{
  result = t(sigEmaxGrad(dose=x, 
                         e0=theta[1], eMax=theta[2], ed50=theta[3], h=theta[4]))
  return(result)
}


aeq_funcsimultaneous7 = function(x, design, gradient, pillarweight, 
                          simultaneouseff, thetapillar, ...){
  # initialize cirterion value
  mysum <- 0
  # calculate equivalence sentence for simultaneous D-optimality (first right site of the equation)
  p <- 4
  right <- p*(sum(pillarweight*simultaneouseff)) 
  
  for(i in 1:7){
    infomat = info(w=design$weights,x=design$supPoints,gradient=gradient,
                   theta= as.numeric(thetapillar[i,]), ...)
    inv_info = solve(infomat)
    current = t(gradient(x,theta= as.numeric(thetapillar[i,]), ...)) %*% inv_info %*% gradient(x, theta= as.numeric(thetapillar[i,]),...)*(pillarweight[i]* simultaneouseff[i])
    mysum=mysum+current
  }
  return(mysum-right)
}



psoOptDesign <- function(crit,control=list(),nPoints=3,Lb=0,Ub=150,xFixed=NULL,wFixed=NULL,xold=NULL,nold=rep(0,length(xold)),n2=rep(0,length(old)),...){
  ###############################################################################
  
  con <- list(numIt = 30, numPart = 300, beta=0.5,gamma=0.7,setRsol=FALSE,setProgressBar=FALSE,OutCritValue=TRUE,intmRes=FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  
  if(!is.null(xFixed)&!is.null(wFixed))
    warning("weights and support points cannot both be predetermined")
  #if weights or support points are predetermined, so is the value of 'nPoints'
  if(!is.null(xFixed)){nPoints<-length(xFixed)}
  if(!is.null(wFixed)){nPoints<-length(wFixed)}
  Blow<-Lb
  Bup<-Ub
  timesteps<-con$numIt
  n<-con$numPart
  betas<-con$beta
  gammas<-con$gamma
  best<-matrix(0,timesteps,nPoints*2)
  
  Lb<-rep(c(0,Lb),each=nPoints)
  Ub<-rep(Ub,2*nPoints)
  
  
  
  #Initializing particles
  
  
  arg<-matrix(0,nPoints*2,n);
  
  if(is.null(wFixed)) {
    arg[1:nPoints,]<-matrix(runif(n*nPoints,0,1),nPoints,n)
    h<-colSums(arg[1:nPoints,])
    
    
    for(i in 1:nPoints){			
      arg[i,]<-arg[i,]/h;
    }
  }
  else{
    arg[1:nPoints,]<-matrix(wFixed,nPoints,n)
  }
  
  if(is.null(xFixed)) {arg[(nPoints+1):(2*nPoints),]<-matrix(runif(n*nPoints,Blow,Bup),nPoints,n)}
  else {arg[(nPoints+1):(2*nPoints),]<-matrix(xFixed,nPoints,n)}
  
  
  zn<-1:n
  
  #==================================================================
  #This is where the actual alogorithm begins
  if(con$setProgressBar) pb<-txtProgressBar(min=0,max=timesteps,style=3)
  for(i in 1:timesteps){
    
    alphas<-gammas^i
    
    for(k in 1:n){
      if(sum(nold)==0) zn[k]<--crit(arg[1:(length(arg[,k])/2),k],arg[(length(arg[,k])/2+1):(length(arg[,k])),k],...)
      else{
        wts<-c(1/(sum(nold)+n2)*nold,(1-sum(nold)/(sum(nold)+n2))*arg[1:(length(arg[,k])/2),k])
        dos<-c(xold,arg[(length(arg[,k])/2+1):(length(arg[,k])),k])
        zn[k]<--crit(wts,dos,...)
      }
    }
    zn_min<-min(zn)
    
    wo<-1:nPoints
    xo<-1:nPoints
    for(k in 1:nPoints){
      wo[k]<-min(arg[k,which.min(zn)])
    }
    for(k in 1:nPoints){
      xo[k]<-min(arg[k+nPoints,which.min(zn)])
    }
    
    zo<-min(zn)
    best[i,]<-c(wo,xo)
    
    #Move the particles
    #the weights...
    if(is.null(wFixed)){
      for(l in 1:n){
        arg[1:nPoints,l]<-arg[1:nPoints,l]*(1-betas)+wo*betas+alphas*rnorm(nPoints,0,0.25)
      }
    }
    # ... and the support points			
    if(is.null(xFixed)){
      for(l in 1:n){
        arg[(nPoints+1):(2*nPoints),l]<-arg[(nPoints+1):(2*nPoints),l]*(1-betas)+xo*betas+alphas*rnorm(nPoints,0,(Bup-Blow)/4)
      }
    }
    
    #Make sure that the particles are still in the search room after having moved
    if(is.null(wFixed))arg[1:nPoints,]<-pmax(arg[1:nPoints,],.1)
    if(is.null(xFixed))arg[(nPoints+1):(2*nPoints),]<-pmax(arg[(nPoints+1):(2*nPoints),],Blow)
    if(is.null(wFixed))arg[1:nPoints,]<-pmin(arg[1:nPoints,],1)
    if(is.null(xFixed))arg[(nPoints+1):(2*nPoints),]<-pmin(arg[(nPoints+1):(2*nPoints),],Bup)
    
    if(is.null(wFixed)){
      h<-colSums(arg[1:nPoints,])
      for(t in 1:nPoints){			
        arg[t,]<-arg[t,]/h;
      }
    }
    if(con$setProgressBar) setTxtProgressBar(pb,i)
    if(con$intmRes)print(best[i,])
  }
  if(con$setProgressBar) close(pb)
  if(con$setRsol){
    print("PSO-Algorithmus has finished. Intermediate results:")
    print(list(weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
    print("start using solnp from package Rsol, please wait...")
    require(Rsolnp, quietly = TRUE)
    eqfun<-function(w) sum(w[1:nPoints])
    fac<-1/crit(wo,xo,...)
    fun<-function(w) -fac*crit(w[1:nPoints],w[(nPoints+1):(2*nPoints)],...)		
    #{ sink("NUL"); z<-round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub,control=list(tol=1e-15))$pars,digits = 5); sink(); }
    z<-round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub)$pars,digits = 5)
    wo<-z[1:nPoints]
    xo<-z[(nPoints+1):(2*nPoints)]
  }
  
  #hin<-function(w) c(w[1:nPoints]+.1,w[(nPoints+1):(2*nPoints)]-Lb,Ub+.1-w[(nPoints+1):(2*nPoints)])
  #constrOptim.nl(c(wo,xo), fun, heq=eqfun, hin=hin)
  #Output the results
  #cat("\n")
  if(con$OutCritValue)
    return(list(value=crit(wo,xo,...),weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
  else return(list(weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
}





psoOptDesign_fixed <-function(crit,control=list(), nPoints=3, Lb=0, Ub= 150, xFixed=NULL,
                       wFixed=NULL, xold=NULL, nold=rep(0,length(xold)),
                       n2=rep(0,length(old)), boundarySup=FALSE, startSup=NULL, 
                       ...){
  ###############################################################################
  
  con <- list(numIt = 30, numPart = 300, beta=0.5, gamma=0.7,
              setRsol=FALSE, setProgressBar=FALSE, 
              OutCritValue=TRUE, intmRes=FALSE)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  
  if(!is.null(xFixed)&!is.null(wFixed))
    warning("weights and support points cannot both be predetermined")
  #if weights or support points are predetermined, so is the value of 'nPoints'
  if(!is.null(xFixed)){nPoints<-length(xFixed)}
  if(!is.null(wFixed)){nPoints<-length(wFixed)}
  
  if(!is.null(startSup)){
    if(length(startSup)!=nPoints)
      warning("length of startvector for the support points need to be equal to nPoints")
  }
  
  # set boundary for design space 
  Blow <- Lb
  Bup <- Ub
  timesteps <- con$numIt
  n <- con$numPart
  betas <- con$beta
  gammas <- con$gamma
  best <- matrix(0,timesteps,nPoints*2)
  
  Lb <- rep(c(0,Lb),each=nPoints)
  Ub <- rep(Ub,2*nPoints)
  
  
  # Initializing particles
  arg <- matrix(0,nPoints*2,n);
  
  if(is.null(wFixed)) {
    arg[1:nPoints,] <- matrix(runif(n*nPoints,0,1),nPoints,n)
    h <- colSums(arg[1:nPoints,])
    
    for(i in 1:nPoints){			
      arg[i,] <- arg[i,]/h;
    }
  }
  else{
    arg[1:nPoints,] <- matrix(wFixed,nPoints,n)
  }
  
  # initialize support points for fixed start or overall fixed support points and
  # random initializing 
  if(is.null(xFixed)) {
    if(is.null(startSup)){
      arg[(nPoints+1):(2*nPoints),] <- matrix(runif(n*nPoints, Blow, Bup), nPoints, n)
    }
    else{arg[(nPoints+1):(2*nPoints),] <- matrix(startSup, nPoints,n)}
  }
  else {arg[(nPoints+1):(2*nPoints),] <- matrix(xFixed,nPoints,n)}
  
  # fix boundary support points if boundarySup is TRUE
  if(boundarySup){
    arg[(nPoints+1),] <- rep(Blow, n) 
    arg[(2*nPoints),] <- rep(Bup, n)}
  else{arg[(nPoints+1):(2*nPoints),] <- arg[(nPoints+1):(2*nPoints),]}
  
  zn <- 1:n
  
  #==================================================================
  #This is where the actual alogorithm begins
  if(con$setProgressBar) pb <- txtProgressBar(min=0, max=timesteps, style=3)
  for(i in 1:timesteps){
    
    alphas<-gammas^i
    
    for(k in 1:n){
      if(sum(nold)==0) zn[k] <- -crit(arg[1:(length(arg[,k])/2),k],arg[(length(arg[,k])/2+1):(length(arg[,k])),k],...)
      else{
        wts<-c(1/(sum(nold)+n2)*nold,(1-sum(nold)/(sum(nold)+n2))*arg[1:(length(arg[,k])/2),k])
        dos<-c(xold,arg[(length(arg[,k])/2+1):(length(arg[,k])),k])
        zn[k] <- -crit(wts,dos,...)
      }
    }
    zn_min <- min(zn)
    
    wo <- 1:nPoints
    xo <- 1:nPoints
    for(k in 1:nPoints){
      wo[k] <- min(arg[k,which.min(zn)])
    }
    for(k in 1:nPoints){
      xo[k] <- min(arg[k+nPoints,which.min(zn)])
    }
    
    zo <- min(zn)
    best[i,] <- c(wo,xo)
    
    #Move the particles
    #the weights...
    if(is.null(wFixed)){
      for(l in 1:n){
        arg[1:nPoints,l] <- arg[1:nPoints,l]*(1-betas)+wo*betas+alphas*rnorm(nPoints,0,0.25)
      }
    }
    # ... and the support points			
    if(is.null(xFixed)){
      if(boundarySup==FALSE){
        for(l in 1:n){
          arg[(nPoints+1):(2*nPoints),l] <- arg[(nPoints+1):(2*nPoints),l]*(1-betas)+
            xo*betas+alphas*rnorm(nPoints,0,(Bup-Blow)/4)
        }
      }
      else{for(l in 1:n){
        arg[(nPoints+2):(2*nPoints-1),l] <- arg[(nPoints+2):(2*nPoints-1),l]*(1-betas)+
          xo[2:(nPoints-1)]*betas+alphas*rnorm(nPoints-2,0,(Bup-Blow)/4)
      }
      }
    }
    
    #Make sure that the particles are still in the search room after having moved
    if(is.null(wFixed))arg[1:nPoints,] <- pmax(arg[1:nPoints,],.1)
    if(is.null(xFixed))arg[(nPoints+1):(2*nPoints),] <- pmax(arg[(nPoints+1):(2*nPoints),],Blow)
    if(is.null(wFixed))arg[1:nPoints,] <- pmin(arg[1:nPoints,],1)
    if(is.null(xFixed))arg[(nPoints+1):(2*nPoints),] <- pmin(arg[(nPoints+1):(2*nPoints),],Bup)
    
    if(is.null(wFixed)){
      h <- colSums(arg[1:nPoints,])
      for(t in 1:nPoints){			
        arg[t,] <- arg[t,]/h;
      }
    }
    if(con$setProgressBar) setTxtProgressBar(pb,i)
    if(con$intmRes)print(best[i,])
  }
  if(con$setProgressBar) close(pb)
  if(con$setRsol){
    print("PSO-Algorithmus has finished. Intermediate results:")
    print(list(weights=round(wo,digits=5),supPoints=round(xo,digits=5)))
    print("start using solnp from package Rsol, please wait...")
    require(Rsolnp, quietly = TRUE)
    eqfun<-function(w) sum(w[1:nPoints])
    fac <- 1/crit(wo,xo,...)
    fun <- function(w) -fac*crit(w[1:nPoints],w[(nPoints+1):(2*nPoints)],...)		
    #{ sink("NUL"); z<-round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub,control=list(tol=1e-15))$pars,digits = 5); sink(); }
    z <- round(solnp(pmin(pmax(c(wo,xo),Blow+1e-15),Bup-1e-15),fun, eqfun, eqB = 1,LB=Lb,UB=Ub)$pars,digits = 5)
    wo <- z[1:nPoints]
    xo <- z[(nPoints+1):(2*nPoints)]
  }
  
  #hin<-function(w) c(w[1:nPoints]+.1,w[(nPoints+1):(2*nPoints)]-Lb,Ub+.1-w[(nPoints+1):(2*nPoints)])
  #constrOptim.nl(c(wo,xo), fun, heq=eqfun, hin=hin)
  #Output the results
  #cat("\n")
  if(con$OutCritValue)
    return(list(value=crit(wo,xo,...), weights=round(wo,digits=5),
                supPoints=round(xo,digits=5)))
  else return(list(weights=round(wo,digits=5),
                   supPoints=round(xo,digits=5)))
}

