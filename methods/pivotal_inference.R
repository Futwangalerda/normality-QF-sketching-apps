
### pivotal inference from sketching estimators of different types,
### including i.i.d., uniform orthogonal, and Hadamard
### for Countsketch, the pivotal result is similar with iidGaussian
pivo_iid_val<-function(k,n,sX,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  lambda<-svd(sX)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmak<-v[k]
  center<-sigmak^2
  est_v<-2*sigmak^4
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}

pivo_iid_projvec<-function(k,n,sX,alpha){
  m<-nrow(sX);p<-ncol(sX)
  lambda<-svd(sX)$d
  RightSingVecs<-svd(sX)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  center<-sqrt(sum(c*RSV[,k])^2)
  ProjVar<-0
  for(i in 1:p){
    if(i!=k){
      ProjVar<-ProjVar+sigma[i]*sigma[k]/(sigma[i]-sigma[k])^2*sum(c*RSV[,i])^2
    }
  }
  est_v<-ProjVar
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}

pivo_haar_val<-function(k,n,sX,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  lambda<-svd(sX)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmak<-v[k]
  center<-sigmak^2
  est_v<-2*(1-gamma)*sigmak^4
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}

pivo_srht_val<-function(k,n,sX,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  lambda<-svd(sX)$d
  v<-sort(lambda,decreasing = TRUE)
  sigmak<-v[k]
  center<-sigmak^2
  est_v<-3*(1-gamma)*sigmak^4
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}

### for projected singular vector, the haar pivotal CI is equal to the hadamard
pivo_srht_projvec<-function(k,n,sX,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  lambda<-svd(sX)$d
  RightSingVecs<-svd(sX)$v
  index<-order(lambda,decreasing = TRUE)
  sigma<-lambda[index]^2
  RSV<-RightSingVecs[,index]
  center<-sqrt(sum(c*RSV[,k])^2)
  ProjVar<-0
  for(i in 1:p){
    if(i!=k){
      ProjVar<-ProjVar+sigma[i]*sigma[k]/(sigma[i]-sigma[k])^2*sum(c*RSV[,i])^2
    }
  }
  est_v<-(1-gamma)*ProjVar
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  
  list(conf=c(center-rb,center-lb))
}