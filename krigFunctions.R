# Copyright KNMI 2021
# MIT License file attached
# jouke.de.baar@knmi.nl

ni = 50
order = 5

# kriging functions
corrMatrix <- function(x1,x2,dem,noisey,sigma,hyper)
{
  H = corrLag(x1,x2,dem,noisey,sigma,hyper)
  A = corrMatrixFast(dem,noisey,sigma,hyper,H$Hhor,H$Hver)
}

corrLag <- function(x1,x2,dem,noisey,sigma,hyper)
{
  N1 = size(x1)
  N2 = size(x2)
  
  # horizontal lag
  X1 = meshgrid(x2[,1,drop=FALSE],x1[,1,drop=FALSE])
  X2 = meshgrid(x2[,2,drop=FALSE],x1[,2,drop=FALSE])
  n1 = size(X1$Y)[1]
  n2 = size(X1$Y)[2]
  lon1 = matrix(X1$Y,n1*n2,1)
  lon2 = matrix(X1$X,n1*n2,1)
  lat1 = matrix(X2$Y,n1*n2,1)
  lat2 = matrix(X2$X,n1*n2,1)
  Hhor = distCosine(cbind(lon1,lat1),cbind(lon2,lat2))
  Hhor = matrix(Hhor,n1,n2)
  
  
  Hver = 0*Hhor
  if(hyper[1]>0){
    for(j in 1:N1[1]){
      for(k in 1:N2[1]){
        loni = linspace(x1[j,1],x2[k,1],ni)
        lati = linspace(x1[j,2],x2[k,2],ni)
        alti = interp2(dem$lat,dem$lon,dem$alt_blur,lati,loni)
        xi = linspace(0,1,ni)
        p = polyfix(xi,alti,order,c(xi[1],xi[ni]),c(alti[1],alti[ni]))
        alti = polyval(p,xi)
        Hver[j,k] = sum(abs(alti[2:ni]-alti[1:(ni-1)]))
      }
    }
  }
  return(list(Hhor=Hhor,Hver=Hver))
}


corrMatrixFast <- function(dem,noisey,sigma,hyper,Hhor,Hver)
{
  H = sqrt(Hhor^2 + (hyper[1]*Hver)^2)
  P = (sigma^2) * ( exp(-0.5*(H^2)/(hyper[2]^2)) )
  if(noisey==0){
    A = P
  } else {
    R = (noisey^2) * eye(size(P)[1])
    A = R + P
  }
  return(A)
}

# kriging: posterior
krigPost <- function(x,y,dem,relnoisey,hyper,xi)
{
  id = is.finite(y)
  x = x[id,,drop=FALSE]
  y = y[id,,drop=FALSE]
  N = length(y)
  sigma = std(y)
  noisey = sigma*relnoisey
  A = corrMatrix(x,x,dem,noisey,sigma,hyper)
  Apd = nearPD(A)
  A = as.matrix(Apd$mat)
  mu = mean(y)
  yn = y - mu
  yn0 = solve(A,yn)
  invA = inv(A)
  ni = size(xi)[1]
  yni = matrix(NA,ni,1)
  ui = matrix(NA,ni,1)
  progress = 1
  for(k in 1:ni){
    if(ni>1000){
      if((100*k/ni)>progress){
        #print(progress)
        progress = progress+1
      }
    }
    b = corrMatrix(xi[k,,drop=FALSE],x,dem,0,sigma,hyper)
    yni[k,] = b %*% yn0
    ui[k,] = sqrt(max(sigma^2 - b%*%invA%*%t(b),0)+noisey^2)
  }
  yi = mu + yni
  return(list(yi=yi,ui=ui))
}

# kriging: LOOCV
krigLOOCVfast <- function(x,y,dem,relnoisey,hyper,Hhor,Hver)
{
  id = is.finite(y)
  x = x[id,,drop=FALSE]
  y = y[id,,drop=FALSE]
  N = length(y)
  sigma = std(y)
  noisey = sigma*relnoisey
  A = corrMatrixFast(dem,noisey,sigma,hyper,Hhor,Hver)
  Apd = nearPD(A)
  A = as.matrix(Apd$mat)
  b = corrMatrixFast(dem,0,sigma,hyper,Hhor,Hver)
  se = matrix(NA,N,1)
  rele = matrix(NA,N,1)
  for(k in 1:N){
    mu = mean(y[-k,,drop=FALSE])
    yn = y[-k,,drop=FALSE]-mu
    Aloc = A[-k,-k,drop=FALSE]
    #Apd = nearPD(Aloc)
    #Aloc = as.matrix(Apd$mat)
    yn0 = solve(Aloc,yn)
    yi = mu + b[k,-k,drop=FALSE]%*%yn0
    covi = sigma^2 - b[k,-k,drop=FALSE]%*%solve(Aloc,t(b[k,-k,drop=FALSE]))
    se[k] = (yi-y[k])^2
    rele[k] = (yi-y[k])/sqrt(max(covi+noisey^2,1e-4))
  }
  rmseyi = sqrt(mean(se))
  nerr = length(rele)
  xcdfObs = sort(rele)
  ycdf = linspace(0.5/nerr,1-0.5/nerr,nerr)
  xcdfOpt = qnorm(ycdf,mean=0,sd=1)
  rmsecdf = mean((xcdfObs-xcdfOpt)^2)
  return(list(rmseyi=rmseyi,rmsecdf=rmsecdf))
}