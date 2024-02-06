### Figure 6
library(mvtnorm)
library(phangorn)
library(PropCIs)
library(ggplot2)
library(tidyr)
library(patchwork)
library(tidyverse)
library(MASS)
setwd('C:/Users/Desktop/test')

### Case 2
set.seed(15)
n=2^12;p=25
genesig<-function(t,n){
  A<-matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      A[i,j]=t^(abs(i-j))
    }
  }
  return(A)
}
sigX<-2*genesig(0.5,p)
W<-rmvt(n, sigX, df = 2)
U<-qr.Q(qr(W));V<-svd(matrix(rnorm(p*p),p))$v
D<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
X<-U%*%D%*%V
beta<-c(rep(1,floor(0.2*p)),0.1*rep(1,p-2*floor(0.2*p)),rep(1,floor(0.2*p)))
y<-X%*%beta+as.vector(rnorm(n,sd=0.01))
squarelength_beta<-sum(beta^2)

sim=2000
grid_m = seq(1600,3600,200)
c<-c(1,rep(0,p-1))
alpha<-0.1

### Hadamard sketching
padding<-function(X,y){
  m<-nrow(X)
  if(ceiling(log(m,2))>log(m,2)){
    m1<-floor(log(m,2))+1
    padX<-rbind(X,matrix(0,2^m1-m,p))
    pady<-append(y,rep(0,2^m1-m))
  }
  else{padX<-X;pady<-y}
  return(list(padX=padX,pady=pady))
}

ols_SRHT<-function(m,c,X,y,partial=0){
  p<-ncol(X)
  pad<-padding(X,y)
  X1<-pad$padX;y1<-pad$pady
  n1<-nrow(X1)
  gamma<-m/n1
  SXy<-apply(sample(c(1,-1),n1,replace=TRUE,prob=c(0.5,0.5))*cbind(X1,y1),2,fhm)[which(rbinom(n1,1,gamma)!=0),]/sqrt(m)
  if(partial==0){g<-qr.solve(SXy[,1:p],SXy[,p+1])}
  else{SX<-SXy[,1:p];M<-solve(t(SX)%*%SX);g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g),r3=SXy))
}

getMLE <- function(x, y, w) {
  return(ginv(x*sqrt(w)) %*% (y*sqrt(w)));
}

####################
##theoretical asymptotic variances for Hadamard sketching, 
#for sketch-and-solve estimators
ep <- y - X %*% beta
XtX_inv<-solve(t(X)%*%X)
th_v<-(1-grid_m/n)*sum(ep^2)*sum((XtX_inv%*%c)*c) 

#for partial sketching estimators
th_v2<-(1-grid_m/n)*(sum((y-ep)^2)*sum((XtX_inv%*%c)*c)+2*sum(c*beta)^2)

#################################
### simulations

t_srht<-matrix(0,sim,length(grid_m))
r_srht<-matrix(0,sim,length(grid_m)*p)
for(j in 1:length(grid_m)){
  for(i in 1:sim){
    start_time<-Sys.time()
    r_srht[i,((j-1)*p+1):(j*p)]<-ols_SRHT(grid_m[j],c,X,y,partial=0)$r1
    end_time<-Sys.time()
    t_srht[i,j]<-as.numeric(difftime(end_time, start_time, units = "secs"))
  }
}

sampledecision <- function(pi){
  pi <- matrix(pi,ncol=1)
  apply(pi,1,rbinom,n=1,size=1)
}

### Optimal Subsampling Algorithm from Yu et al. (2022)
AlgTwoStp <- function(r1=r1, r2=r2) {
  if (r2 == 0) {
    idx <- 1:n
    pi <- rep((r1/n),n)
    decision <- rbinom(n,rep(1,n),prob=pi)
    idx.simp <- idx[decision==1]
    x.simp <- X[idx.simp,]
    y.simp <- y[idx.simp]
    fit.simp <- ginv(x.simp) %*% y.simp
    beta.simp <- fit.simp
    return(simp=beta.simp)
  }
  if (r2 != 0) {
    idx <- 1:n
    idx.prop <- sample(1:n, r1, T)
    x.prop <- X[idx.prop,]
    y.prop <- y[idx.prop]
    pinv.prop <- rep(n,r1)
    fit.prop <- getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop
    psi.dot  <-  X %*% beta.prop
    ## mVc
    PI.mVc <- abs(y - psi.dot) * rowSums(X^2)
    PI.mVc <- PI.mVc / sum(PI.mVc)
    PI.mVc1 <- abs((1-0.01)*PI.mVc + 0.01/n+1/(r2+1e-6))/2-abs((1-0.01)*PI.mVc + 0.01/n-1/(r2+1e-6))/2
    decision <- rbinom(n,rep(1,n),prob=r2*PI.mVc1)
    idx.mVc <- idx[decision==1]
    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- y[c(idx.mVc, idx.prop)]
    pinv.mVc <- c(1 / PI.mVc1[idx.mVc], pinv.prop)
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    return(fit.mVc)
  }
}

### Different Optimal Subsampling Methods
r1 <- 800
r_optsubsamp_800<-matrix(0,sim,length(grid_m)*p)
t_optsubsamp_800<-matrix(0,sim,length(grid_m))
for(j in 1:length(grid_m)){
  for(i in 1:sim){
    start_time<-Sys.time()
    r2 <- grid_m[j]-r1
    r_optsubsamp_800[i,((j-1)*p+1):(j*p)]<-AlgTwoStp(r1, r2)
    end_time<-Sys.time()
    t_optsubsamp_800[i,j]<-as.numeric(difftime(end_time, start_time, units = "secs"))
  }
}
r1 <- 1200
r_optsubsamp_1200<-matrix(0,sim,length(grid_m)*p)
t_optsubsamp_1200<-matrix(0,sim,length(grid_m))
for(j in 1:length(grid_m)){
  for(i in 1:sim){
    start_time<-Sys.time()
    r2 <- grid_m[j]-r1
    r_optsubsamp_1200[i,((j-1)*p+1):(j*p)]<-AlgTwoStp(r1, r2)
    end_time<-Sys.time()
    t_optsubsamp_1200[i,j]<-as.numeric(difftime(end_time, start_time, units = "secs"))
  }
}
r1 <- 1400
r_optsubsamp_1400<-matrix(0,sim,length(grid_m)*p)
t_optsubsamp_1400<-matrix(0,sim,length(grid_m))
for(j in 1:length(grid_m)){
  for(i in 1:sim){
    start_time<-Sys.time()
    r2 <- grid_m[j]-r1
    r_optsubsamp_1400[i,((j-1)*p+1):(j*p)]<-AlgTwoStp(r1, r2)
    end_time<-Sys.time()
    t_optsubsamp_1400[i,j]<-as.numeric(difftime(end_time, start_time, units = "secs"))
  }
}

### Uniform Subsampling
r_unifsubsamp<-matrix(0,sim,length(grid_m)*p)
t_unifsubsamp<-matrix(0,sim,length(grid_m))
for(j in 1:length(grid_m)){
  for(i in 1:sim){
    start_time<-Sys.time()
    fit.alg <- AlgTwoStp(grid_m[j], 0)
    r_unifsubsamp[i,((j-1)*p+1):(j*p)]<-fit.alg
    end_time<-Sys.time()
    t_unifsubsamp[i,j]<-as.numeric(difftime(end_time, start_time, units = "secs"))
  }
}

#################################################
### display the MSE

srht_var<-apply(r_srht, 2, var)
srht_mse<-rep(0, length(grid_m))
for (j in 1:length(grid_m)){
  srht_mse[j]<-(grid_m[j]/p)*sum(srht_var[((j-1)*p+1):(j*p)])/squarelength_beta
}
optsubsamp_1200_var<-apply(r_optsubsamp_1200, 2, var)
optsubsamp_1200_mse<-rep(0, length(grid_m))
for (j in 1:length(grid_m)){
  optsubsamp_1200_mse[j]<-(grid_m[j]/p)*sum(optsubsamp_1200_var[((j-1)*p+1):(j*p)])/squarelength_beta
}
optsubsamp_1400_var<-apply(r_optsubsamp_1400, 2, var)
optsubsamp_1400_mse<-rep(0, length(grid_m))
for (j in 1:length(grid_m)){
  optsubsamp_1400_mse[j]<-(grid_m[j]/p)*sum(optsubsamp_1400_var[((j-1)*p+1):(j*p)])/squarelength_beta
}
optsubsamp_800_var<-apply(r_optsubsamp_800, 2, var)
optsubsamp_800_mse<-rep(0, length(grid_m))
for (j in 1:length(grid_m)){
  optsubsamp_800_mse[j]<-(grid_m[j]/p)*sum(optsubsamp_800_var[((j-1)*p+1):(j*p)])/squarelength_beta
}
unifsubsamp_var<-apply(r_unifsubsamp, 2, var)
unifsubsamp_mse<-rep(0, length(grid_m))
for (j in 1:length(grid_m)){
  unifsubsamp_mse[j]<-(grid_m[j]/p)*sum(unifsubsamp_var[((j-1)*p+1):(j*p)])/squarelength_beta
}

t_srht_mean<-apply(t_srht, 2, mean)
t_srht_var<-apply(t_srht, 2, var)
t_optsubsamp_800_mean<-apply(t_optsubsamp_800, 2, mean)
t_optsubsamp_800_var<-apply(t_optsubsamp_800, 2, var)
t_optsubsamp_1200_mean<-apply(t_optsubsamp_1200, 2, mean)
t_optsubsamp_1200_var<-apply(t_optsubsamp_1200, 2, var)
t_optsubsamp_1400_mean<-apply(t_optsubsamp_1400, 2, mean)
t_optsubsamp_1400_var<-apply(t_optsubsamp_1400, 2, var)
t_unifsubsamp_mean<-apply(t_unifsubsamp, 2, mean)
t_unifsubsamp_var<-apply(t_unifsubsamp, 2, var)

Methods<-c(rep("hadamard",length(grid_m)),
                   rep("optsubsamp_800",length(grid_m)),
                   rep("optsubsamp_1200",length(grid_m)),
                   rep("optsubsamp_1400",length(grid_m)),
                   rep("unifsubsamp",length(grid_m)))
standard_mse<-c((srht_mse), (optsubsamp_800_mse), (optsubsamp_1200_mse), (optsubsamp_1400_mse), (unifsubsamp_mse))
time_mean<-c(t_srht_mean, t_optsubsamp_800_mean,t_optsubsamp_1200_mean,t_optsubsamp_1400_mean,t_unifsubsamp_mean)
time_var<-c(t_srht_var, t_optsubsamp_800_var,t_optsubsamp_1200_var,t_optsubsamp_1400_var,t_unifsubsamp_var)
sum_grid_m<-rep(grid_m, 5)
data_sum<-data.frame(sum_grid_m, time_mean, time_var, standard_mse, Methods)

ggplot(data_sum, aes(time_mean, standard_mse, color=Methods, shape=Methods))+
  geom_point(size = 1.5)+
  geom_errorbar(aes(ymin = standard_mse - qnorm(1-alpha/2,sd=sqrt(2*standard_mse^2/(sum_grid_m-1))), ymax = standard_mse - qnorm(alpha/2,sd=sqrt(2*standard_mse^2/(sum_grid_m-1)))), width = 0.0001, alpha=0.2)+
  geom_errorbarh(aes(xmin= time_mean - qnorm(1-alpha/2,sd=sqrt(time_var/sim)), xmax=time_mean - qnorm(alpha/2,sd=sqrt(time_var/sim))), height=0.05, alpha=0.2)+
  theme_gray()+
    theme(
      legend.position = "right", legend.key.size = unit(0.8, 'cm'), legend.text=element_text(size=12),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(color = "grey20", size = 15),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(color = "grey20", size = 15)
    ) +
    xlab('Time') +
    ylab('Standard MSE')+
    ggtitle('Complete Sketching, Case 2')

file_path<-"C:/Users/Desktop/test/"
ggsave("n4096p25Case2_var_time_comp_wallclock.png", path=file_path, width = 6, height = 4, units = "in", dpi = 300)

