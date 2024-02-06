library(phangorn)

#iid sketching
Estising_iid<-function(i,m,X,c,type=2){
  n<-nrow(X)
  p<-ncol(X)
  if(type==1){S<-matrix(rnorm(m*n),m,n)/sqrt(m)}
  if(type==2){nbinom<-6
  p0<-0.5+1/(2*sqrt(3))
  S<-matrix(rbinom(m*n,nbinom,p0)-nbinom*p0,m,n)/sqrt(m)}
  P<-S%*%X
  lambda<-svd(P)$d
  RightSingVecs<-svd(P)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  return(list(lambdai=sigma[i], squareprojveci=abs(sum(c*RSV[,i])), matrix=P))
}
#random orthogonal sketching
Generate_Haar<-function(m,n){
  O<-matrix(rnorm(m*n),m,n)
  S<-t(svd(O)$v)*sqrt(n/m)
  return(S)
}
Estising_Haar<-function(i,m,X,c){
  n<-nrow(X)
  p<-ncol(X)
  S_haar<-Generate_Haar(m,n)
  P<-S_haar%*%X
  lambda<-svd(P)$d
  RightSingVecs<-svd(P)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  return(list(lambdai=sigma[i], squareprojveci=abs(sum(c*RSV[,i])), matrix=P))
}

### Hadamard sketching
padding<-function(X){
  m<-nrow(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
  }
  else{padX<-X}
  return(padX)
}

Estising_SRHT<-function(i,m,X,c){
  n<-nrow(X)
  p<-ncol(X)
  X1<-padding(X)
  n1<-nrow(X1)
  gamma<-m/n1
  SX<-apply(sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))*X1,2,fhm)[which(rbinom(n1,1,gamma)!=0),]/sqrt(m)
  lambda<-svd(SX)$d
  RightSingVecs<-svd(SX)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  return(list(lambdai=sigma[i], squareprojveci=abs(sum(c*RSV[,i])), matrix=SX))
}
### CountSketch
Estising_CountSketch<-function(i,m,X,c){
  n<-nrow(X)
  p<-ncol(X)
  SX<-matrix(0,m,p)
  randSigns<-sample(c(1,-1),n,replace=TRUE,prob=c(0.5,0.5))
  hash<-sample(m,n,replace=TRUE)
  for(j in 1:n){
    a <- X[j, ]
    h <- hash[j]
    g <- randSigns[j]
    SX[h, ] <- SX[h, ] + g * a
  }
  lambda<-svd(SX)$d
  RightSingVecs<-svd(SX)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  return(list(lambdai=sigma[i], squareprojveci=abs(sum(c*RSV[,i])), matrix=SX))
}
### Uniform Subsampling Method
Estising_UnifSubsamp<-function(i,m,X,c){
  n<-nrow(X)
  gamma<-m/n
  SX <- X[which(rbinom(n,1,gamma)!=0), ]
  lambda<-svd(SX)$d
  RightSingVecs<-svd(SX)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  return(list(lambdai=sigma[i], squareprojveci=abs(sum(c*RSV[,i])), matrix=SX))
}