### Figure 4
### Display the variance of different sketching methods

library(mvtnorm)
library(phangorn)
library(PropCIs)
library(ggplot2)
library(tidyr)
library(patchwork)
library(nycflights13)
library(tidyverse)
library(MASS)
setwd('C:/Users/ldwan/Desktop/test111111')

### NYCflight Dataset
flt<-flights %>% na.omit() #%>% filter(10<= month & month<=12)
head(flt)
A <- merge(flt, weather)
B <- merge(A, planes, by = "tailnum")
C <- merge(B, airlines, by = "carrier")
variate <- c("arr_delay", "month", "day", "dep_time", "sched_dep_time", "dep_delay", "arr_time",
             "sched_arr_time", "air_time", "distance", "hour", "temp", "dewp", "humid",
             "wind_dir", "wind_speed", "wind_gust", "precip", "pressure", "visib", "year.y", "seats")
dat <- C[, variate]
dat <- na.omit(dat)
X1<-as.matrix(dat[,2:22])
y1<-as.matrix(dat[,1])
X<-scale(X1)
y<-scale(y1)

### Case 1
# set.seed(15)
# n=2^11;p=15
# d<-1/(1:p)
# D<-diag(d)
# O1<-matrix(rnorm(n*p),n,p)
# W<-svd(O1)$u
# O2<-matrix(rnorm(p*p),p,p)
# U<-qr.Q(qr(O2))
# X<-W %*% D %*% t(U)
# y<-runif(n,0,1)
# set.seed(NULL)

### Case 2
# genesig<-function(t,n){
#   A<-matrix(0,n,n)
#   for(i in 1:n){
#     for(j in 1:n){
#       A[i,j]=t^(abs(i-j))
#     }
#   }
#   return(A)
# }
# sigX<-2*genesig(0.5,p)
# W<-rmvt(n, sigX, df = 2)
# U<-qr.Q(qr(W));V<-svd(matrix(rnorm(p*p),p))$v
# D<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
# X<-U%*%D%*%V
# beta<-c(rep(1,floor(0.2*p)),0.1*rep(1,p-2*floor(0.2*p)),rep(1,floor(0.2*p)))
# y<-X%*%beta+as.vector(rnorm(n,sd=0.01))

### Case 3
# X1<-matrix(rnorm(n/2*p),n/2,p)
# X2<-matrix(rnorm(n/2*p,5,1),n/2,p)
# X<-rbind(X1,X2)
# y<-runif(n,0,1)

Xty<-t(X)%*%y
n<-nrow(X)
p<-ncol(X)
beta<-qr.solve(X,y)

sim=500
grid_m = seq(2000,10000,1000)
c<-c(1,rep(0,p-1))
alpha<-0.1

sampleDist<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/3,1/3,1/3))*sqrt(3/2)
}

sampleDist2<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/6,2/3,1/6))*sqrt(3)
}

sampleDistsparse<-function(n){
  sample(x=c(-1,0,1),n,replace=T,prob=c(1/20,9/10,1/20))*sqrt(10)
}

# different types of S have different kurtosis
ols_iid<-function(m,c,X,y,type=4,partial=0){
  n<-nrow(X)
  if(type==1){S<-matrix(sampleDist(m*n),m)/sqrt(m)}  
  if(type==2){S<-matrix(runif(m*n,-sqrt(3),sqrt(3)),m)/sqrt(m)}
  if(type==3){S<-matrix(sampleDist2(m*n),m)/sqrt(m)}
  if(type==4){S<-matrix(rnorm(m*n),m)/sqrt(m)}
  if(type==5){S<-matrix(rt(m*n,10),m,n)/sqrt(10/8)/sqrt(m)}
  if(type==6){S<-matrix(sampleDistsparse(m*n),m)/sqrt(m)}
  
  SX<-S%*%X;Sy<-S%*%y
  M<-solve(t(SX)%*%SX)
  if(partial==0){w<-t(SX)%*%Sy;g<-M%*%w}
  else{g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g),r3=cbind(SX,Sy)))
}



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


ols_CountSketch<-function(m,c,X,y,partial=0){
  n <- nrow(X)
  p <- ncol(X)
  Xy<-cbind(X,y)
  SXy <- matrix(0, m, p+1)
  hashedIndices <- sample(m, n, replace = TRUE)
  randSigns <- sample(c(-1, 1), n, replace = TRUE)
  for (j in 1:n) {
    a <- Xy[j, ]
    h <- hashedIndices[j]
    g <- randSigns[j]
    SXy[h, ] <- SXy[h, ] + g * a
  }
  if(partial==0){g<-qr.solve(SXy[,1:p],SXy[,p+1])}
  else{SX<-SXy[,1:p];M<-solve(t(SX)%*%SX);g<-M%*%(t(X)%*%y)}
  return(list(r1=g,r2=sum(c*g),r3=SXy))
}

####################
##theoretical asymptotic variances for Hadamard sketching, 
#for sketch-and-solve estimators
ep <- y - X %*% beta
XtX_inv<-solve(t(X)%*%X)
th_v<-(1-grid_m/n)*sum(ep^2)*sum((XtX_inv%*%c)*c) 

#for partial sketching estimators
th_v2<-(1-grid_m/n)*(sum((y-ep)^2)*sum((XtX_inv%*%c)*c)+2*sum(c*beta)^2)

####################
###theoretical asymptotic variances for i.i.d. sketching estimators
kappa4 = 3
ep<- as.vector(y - X%*%beta)
Xinv<-solve(t(X)%*%X)%*%t(X)
l<-as.vector(t(c)%*%Xinv)
th_v_iid<-(kappa4 - 3)*sum((l*ep)^2)+ sum(l^2)*sum(ep^2)
th_v2_iid<-(kappa4 - 3)*sum((l*(y-ep))^2)+
  (sum(l^2)*sum((y-ep)^2)+sum(c*beta)^2)
#################################
### simulations
# sketch-and-solve for i.i.d., Hadamard, and Countsketch

sim=500
r_iid<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_iid[i,j]<-ols_iid(grid_m[j],c,X,y,type=4)$r2}
}

r_srht<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_srht[i,j]<-ols_SRHT(grid_m[j],c,X,y,partial=0)$r2}
}

r_countsketch<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_countsketch[i,j]<-ols_CountSketch(grid_m[j],c,X,y,partial=0)$r2}
}

r_iid_p<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_iid_p[i,j]<-ols_iid(grid_m[j],c,X,y,type=4,partial=1)$r2}
}



r_srht_p<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_srht_p[i,j]<-ols_SRHT(grid_m[j],c,X,y,partial=1)$r2}
}


r_countsketch_p<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    r_countsketch_p[i,j]<-ols_CountSketch(grid_m[j],c,X,y,partial=1)$r2}
}

#################################################
### display the variances
sketch_size<-as.vector(t(matrix(rep(grid_m,sim),length(grid_m))))
df_var1<-data.frame(sketch_size, iid=as.vector(r_iid)*sqrt(sketch_size),hadamard=as.vector(r_srht)*sqrt(sketch_size),countsketch=as.vector(r_countsketch)*sqrt(sketch_size))
df_var2<-data.frame(sketch_size, iid=as.vector(r_iid_p)*sqrt(sketch_size),hadamard=as.vector(r_srht_p)*sqrt(sketch_size),countsketch=as.vector(r_countsketch_p)*sqrt(sketch_size))

v1<-aggregate(.~sketch_size,df_var1, FUN = stats::var)
v2<-aggregate(.~sketch_size,df_var2, FUN = stats::var)

v1$hadamard_theory<-th_v;v2$hadamard_theory<-th_v2
v1$iid_theory<-th_v_iid;v2$iid_theory<-th_v2_iid
res_vr1<-gather(v1,type,var, -sketch_size)
res_vr2<-gather(v2,type,var, -sketch_size)
res_vr1$var<-log(res_vr1$var,10)
res_vr2$var<-log(res_vr2$var,10)
data_c<-data.frame(x = v1$sketch_size, 
                   y1=log(v1$iid), 
                   y2=log(v1$hadamard), y3=log(v1$countsketch), y4=log(v1$iid_cs_theory), y5=log(v1$hadamard_theory))
data_p<-data.frame(x = v2$sketch_size, 
                   y1=log(v2$iid), 
                   y2=log(v2$hadamard), y3=log(v2$countsketch), y4=log(v2$iid_cs_theory), y5=log(v2$hadamard_theory))


p1f<-ggplot(data_c,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 1.5) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.5) +
  geom_point(aes(y = y2, color = "hadamard",shape="hadamard"), size = 1.5) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.5) +
  geom_point(aes(y = y3, color = "countsketch",shape="countsketch"), size = 1.5) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.5) +
  geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 1.5) +
  geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.5) +
  geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 1.5) +
  geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.5) +
  geom_errorbar(aes(ymin= log(exp(y1) - qnorm(1-alpha/2,sd=sqrt(2*exp(y1)^2/(x-1)))), ymax=log(exp(y1) - qnorm(alpha/2,sd=sqrt(2*exp(y1)^2/(x-1)))), color="iid"), width=50)+
  geom_errorbar(aes(ymin= log(exp(y2) - qnorm(1-alpha/2,sd=sqrt(2*exp(y2)^2/(x-1)))), ymax=log(exp(y2) - qnorm(alpha/2,sd=sqrt(2*exp(y2)^2/(x-1)))), color="hadamard"), width=50)+
  geom_errorbar(aes(ymin= log(exp(y3) - qnorm(1-alpha/2,sd=sqrt(2*exp(y3)^2/(x-1)))), ymax=log(exp(y3) - qnorm(alpha/2,sd=sqrt(2*exp(y3)^2/(x-1)))), color="countsketch"), width=50)+
  scale_shape_manual(values = c("iid" = 1, "hadamard" = 2, "countsketch" = 3, "iid_cs_theory" = 4,"hadamard_theory"=5)) +
  scale_color_manual(values = c("iid" = 'green', "hadamard" = 'red', "countsketch" = 'blue', "iid_cs_theory" = 'black',"hadamard_theory"='#FFBB22')) +
  scale_linetype_manual(values=c("iid" = 'solid', "hadamard" = 'longdash', "countsketch" = 'twodash', "iid_cs_theory" = 'dashed',"hadamard_theory"='dotted'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_gray()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(color = "grey20", size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(color = "grey20", size = 20)
  ) +
  xlab('m') +
  ylab('Log Variance')+
  ggtitle("Complete Sketching")

p1p<-ggplot(data_p,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 1.5) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.5) +
  geom_point(aes(y = y2, color = "hadamard",shape="hadamard"), size = 1.5) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.5) +
  geom_point(aes(y = y3, color = "countsketch",shape="countsketch"), size = 1.5) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.5) +
  geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 1.5) +
  geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.5) +
  geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 1.5) +
  geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.5) +
  geom_errorbar(aes(ymin= log(exp(y1) - qnorm(1-alpha/2,sd=sqrt(2*exp(y1)^2/(x-1)))), ymax=log(exp(y1) - qnorm(alpha/2,sd=sqrt(2*exp(y1)^2/(x-1)))), color="iid"), width=50)+
  geom_errorbar(aes(ymin= log(exp(y2) - qnorm(1-alpha/2,sd=sqrt(2*exp(y2)^2/(x-1)))), ymax=log(exp(y2) - qnorm(alpha/2,sd=sqrt(2*exp(y2)^2/(x-1)))), color="hadamard"), width=50)+
  geom_errorbar(aes(ymin= log(exp(y3) - qnorm(1-alpha/2,sd=sqrt(2*exp(y3)^2/(x-1)))), ymax=log(exp(y3) - qnorm(alpha/2,sd=sqrt(2*exp(y3)^2/(x-1)))), color="countsketch"), width=50)+
  scale_shape_manual(values = c("iid" = 1, "hadamard" = 2, "countsketch" = 3, "iid_cs_theory" = 4,"hadamard_theory"=5)) +
  scale_color_manual(values = c("iid" = 'green', "hadamard" = 'red', "countsketch" = 'blue', "iid_cs_theory" = 'black',"hadamard_theory"='#FFBB22')) +
  scale_linetype_manual(values=c("iid" = 'solid', "hadamard" = 'longdash', "countsketch" = 'twodash', "iid_cs_theory" = 'dashed',"hadamard_theory"='dotted'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_gray()+
  theme(
    legend.position = "none",
    plot.title = element_text(size = 22, hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(color = "grey20", size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(color = "grey20", size = 20)
  ) +
  xlab('m') +
  ylab('Log Variance')+
  ggtitle("Partial Sketching")

