### Figure 4
### Compare the variances of different sketching ls solutions
library(mvtnorm)
library(phangorn)
library(ggplot2)
library(tidyr)
library(patchwork)
setwd('C:/Users/Desktop/test/')


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


n=2^11;p=15
grid_m=seq(200,1600,200)
c<-c(1,rep(0,p-1))
alpha<-0.1
###############################################################
### Case 1 
set.seed(15)
n=2^11;p=15
d<-1/(1:p)
D<-diag(d)
O1<-matrix(rnorm(n*p),n,p)
W<-svd(O1)$u
O2<-matrix(rnorm(p*p),p,p)
U<-qr.Q(qr(O2))
X<-W %*% D %*% t(U)
y<-runif(n,0,1)
set.seed(NULL)
ls<-qr.solve(X,y)
sum((X%*%ls)^2)/sum(y^2)
ls[1]


### Case 2
# set.seed(15)
# n=2^11;p=15
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
# ls<-qr.solve(X,y)
# ls[1]

### Case 3
# set.seed(15)
# n=2^11;p=15
# X1<-matrix(rnorm(n/2*p),n/2,p)
# X2<-matrix(rnorm(n/2*p,5,1),n/2,p)
# X<-rbind(X1,X2)
# #y<-1/(1:n)
# y<-runif(n,0,1)
# ls<-qr.solve(X,y)
# ls[1]

#signal-to-noise ratio
lsbeta<-qr.solve(X,y)
MSS<-sum((X%*%lsbeta)^2)
TSS<-sum(y^2)
snr<-MSS/(TSS-MSS);snr

####################
##theoretical asymptotic variances for Hadamard sketching, 
#for sketch-and-solve estimators
ep<-y-X%*%lsbeta
XtX_inv<-solve(t(X)%*%X)
th_v<-(1-grid_m/n)*sum(ep^2)*sum((XtX_inv%*%c)*c) 

#for partial sketching estimators
th_v2<-(1-grid_m/n)*(sum((y-ep)^2)*sum((XtX_inv%*%c)*c)+2*sum(c*lsbeta)^2)

####################
###theoretical asymptotic variances for i.i.d. sketching estimators
kappa4 = 3
ep<- as.vector(y - X%*%lsbeta)
Xinv<-solve(t(X)%*%X)%*%t(X)
l<-as.vector(t(c)%*%Xinv)
th_v_iid<-(kappa4 - 3)*sum((l*ep)^2)+ sum(l^2)*sum(ep^2)
th_v2_iid<-(kappa4 - 3)*sum((l*(y-ep))^2)+
  (sum(l^2)*sum((y-ep)^2)+sum(c*lsbeta)^2)


#################################
### simulations
# sketch-and-solve for i.i.d., Hadamard, and Countsketch

### Complete Sketching
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


### Partial Sketching

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
v1$iid_cs_theory<-th_v_iid;v2$iid_cs_theory<-th_v2_iid
data_c<-data.frame(x = v1$sketch_size, y1=log(v1$iid), y2=log(v1$hadamard), y3=log(v1$countsketch), y4=log(v1$iid_cs_theory), y5=log(v1$hadamard_theory))
data_p<-data.frame(x = v2$sketch_size, y1=log(v2$iid), y2=log(v2$hadamard), y3=log(v2$countsketch), y4=log(v2$iid_cs_theory), y5=log(v2$hadamard_theory))
write.csv(data_c,"C:/Users/ldwan/Desktop/test111111/n2048p15_ols_Case1_complete_var.csv")
write.csv(data_p,"C:/Users/ldwan/Desktop/test111111/n2048p15_ols_Case1_partial_var.csv")


p1f<-ggplot(data_c,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape= 'iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "hadamard",shape= "hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "countsketch",shape= "countsketch"), size = 2) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
  geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 2) +
  geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.2) +
  geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 2) +
  geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.2) +
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
    legend.position = "right", legend.key.size = unit(0.8, 'cm'), legend.text=element_text(size=12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(color = "grey20", size = 15)
  ) +
  xlab('m') +
  ylab('Log Variance')+
  #ylab('Coverage')+
  ggtitle('Complete Sketching')

p1p<-ggplot(data_p,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape= 'iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "hadamard",shape= "hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "countsketch",shape= "countsketch"), size = 2) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
  geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 2) +
  geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.2) +
  geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 2) +
  geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.2) +
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
    legend.position = "right", legend.key.size = unit(0.8, 'cm'), legend.text=element_text(size=12),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(color = "grey20", size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(color = "grey20", size = 15)
  ) +
  xlab('m') +
  ylab('Log Variance')+
  ggtitle('Partial Sketching')
