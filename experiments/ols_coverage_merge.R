### Figure 3

library(mvtnorm)
library(phangorn)
library(PropCIs)
library(ggplot2)
library(tidyr)
library(patchwork)
library(nycflights13)

##############################################################
###  Generate sketched data
### i.i.d. sketching
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

### Countsketch
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



#########################################################
### find confidence interval


pivo_iid<-function(c,n,sX,sy,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  gamma<-m/n
  invsX<-solve((t(sX)%*%sX))
  if(partial==0){
    sX_pinv<-invsX%*%t(sX)
    lssk<-sX_pinv%*%sy
    center<-sum(c*lssk)
    sep<-sy-sX%*%lssk
    a<-t(sX_pinv)%*%c
    est_v<-sum((a*sep)^2)
  }
  else{
    lssk<-invsX%*%Xty
    center<-sum(c*lssk)
    sX_pinv<-invsX%*%t(sX)
    sXbeta<-sX%*%lssk
    a<-t(sX_pinv)%*%c
    est_v<-sum((a*sXbeta)^2)
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v))
  lb=qnorm(alpha/2,sd=sqrt(est_v))
  
  list(conf=c(center-rb,center-lb))
}



pivo_hadamard<-function(c,n,sX,sy,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  invsX<-solve((t(sX)%*%sX))
  if(partial==0){
    lssk<-invsX%*%(t(sX)%*%sy)
    center<-sum(c*lssk)
    sep<-sy-sX%*%lssk
    est_v<-(1-gamma)*sum(c*(invsX%*%c))*sum(sep^2)
  }
  else{
    lssk<-invsX%*%Xty
    center<-sum(c*lssk)### for partial case, we assume the knowledge of $X^\top y$
    est_v<-(1-gamma)*((sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+2*(sum(c*lssk))^2))  
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(center-rb,center-lb))
}


pivo_countsketch<-function(c,n,sX,sy,partial=0,alpha){
  m<-nrow(sX);p<-ncol(sX)
  xi<-p/n;gamma<-m/n
  invsX<-solve((t(sX)%*%sX))
  if(partial==0){
    lssk<-invsX%*%(t(sX)%*%sy)
    center<-sum(c*lssk)
    sep<-sy-sX%*%lssk
    est_v<-sum(c*(invsX%*%c))*sum(sep^2)
  }
  else{
    lssk<-invsX%*%Xty
    center<-sum(c*lssk)### for partial case, we assume the knowledge of $X^\top y$
    est_v<-sum((sX%*%lssk)^2)*sum(c*(invsX%*%c))+(sum(c*lssk))^2  
    ### the asymptotic variance of Hadamard partial sketching estimators are slightly different from that of Haar case
  }
  rb=qnorm(1-alpha/2,sd=sqrt(est_v/m))
  lb=qnorm(alpha/2,sd=sqrt(est_v/m))
  list(conf=c(center-rb,center-lb))
}


################################################################
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

### NYCflight Dataset
# flt<-flights %>% na.omit() #%>% filter(10<= month & month<=12)
# head(flt)
# A <- merge(flt, weather)
# B <- merge(A, planes, by = "tailnum")
# C <- merge(B, airlines, by = "carrier")
# variate <- c("arr_delay", "month", "day", "dep_time", "sched_dep_time", "dep_delay", "arr_time",
#              "sched_arr_time", "air_time", "distance", "hour", "temp", "dewp", "humid",
#              "wind_dir", "wind_speed", "wind_gust", "precip", "pressure", "visib", "year.y", "seats")
# dat <- C[, variate]
# dat <- na.omit(dat)
# X1<-as.matrix(dat[,2:22])
# y1<-as.matrix(dat[,1])
# X<-scale(X1)
# y<-scale(y1)

ls<-qr.solve(X,y)
sum((X%*%ls)^2)/sum(y^2)
ls[1]

c<-c(1,rep(0,ncol(X)-1))
Xty<-t(X)%*%y
grid_m<-seq(200,1600,200)

sim=500
pivo_conf_iid_s=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_iid(grid_m[j],c,X,y,type=4,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    pivo_conf_iid_s[i,(2*j-1):(2*j)]<-pivo_iid(c,n,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_iid_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_iid(grid_m[j],c,X,y,type=4,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_iid_pa[i,(2*j-1):(2*j)]<-pivo_iid(c,n,SX,Sy,partial=1,alpha=0.05)$conf
  }
}


Xy<-padding(X,y)
X0<-Xy$padX;y0<-Xy$pady
n0=nrow(X0)
pivo_conf_hadamard_s=matrix(0,sim,2*length(grid_m))

for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_SRHT(grid_m[j],c,X0,y0,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    pivo_conf_hadamard_s[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n0,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_hadamard_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_SRHT(grid_m[j],c,X0,y0,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    lssk<-res$r1
    pivo_conf_hadamard_pa[i,(2*j-1):(2*j)]<-pivo_hadamard(c,n0,SX,Sy,partial=1,alpha=0.05)$conf
  }
}

pivo_conf_countsketch_s=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_CountSketch(grid_m[j],c,X,y,partial=0)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    pivo_conf_countsketch_s[i,(2*j-1):(2*j)]<-pivo_countsketch(c,n,SX,Sy,partial=0,alpha=0.05)$conf
  }
}

pivo_conf_countsketch_pa=matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    res<-ols_CountSketch(grid_m[j],c,X,y,partial=1)
    SX<-res$r3[,1:p];Sy<-res$r3[,p+1]
    pivo_conf_countsketch_pa[i,(2*j-1):(2*j)]<-pivo_countsketch(c,n0,SX,Sy,partial=1,alpha=0.05)$conf
  }
}

conf<-cbind(pivo_conf_iid_s,pivo_conf_hadamard_s,pivo_conf_countsketch_s)
conf_pa<-cbind(pivo_conf_iid_pa,pivo_conf_hadamard_pa,pivo_conf_countsketch_pa)

## Coverage of 90\% intervals for the first coordinate of $\beta_n$,  and 95\% Clopper-Pearson interval for the coverage 
even<-seq(2,ncol(conf),2);odd<-seq(1,ncol(conf)-1,2)
right<-conf[,even];left<-conf[,odd]
accept<-(left<ls[1])&(ls[1]<right)
cov0<-apply(accept,2,sum)
cov<-matrix(cov0,ncol(conf)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio<-cov/sim  #  rows represent different methods, and columns represent different sketching size m
ratio
ci<-function(x){exactci(x,sim,0.95)$conf.int[1:2]} 
int=matrix(0,ncol(conf)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    int[i,(2*j-1):(2*j)]=round(ci(cov[i,j]),digits=3)
  }
}

lower<-int[,seq(1,2*length(grid_m),2)]
upper<-int[,seq(2,2*length(grid_m),2)]


x <- grid_m
df_full <- data.frame(x = grid_m, y1 =ratio[1,], y2=ratio[2,],y3=ratio[3,], y_upper1 =upper[1,] ,
                      y_lower1 = lower[1,],y_upper2=upper[2,],y_lower2=lower[2,],y_upper3=upper[3,],y_lower3=lower[3,])


p1f_s<-ggplot(df_full, aes(x = x)) +
  ylim(0.9, 0.98) +
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 2) +
  geom_point(aes(y = y2, color = "Hadamard",shape="Hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "Hadamard",linetype = 'Hadamard'), linewidth = 2) +
  geom_point(aes(y = y3, color = "CountSketch",shape="CountSketch"), size = 2) +
  geom_line(aes(y = y3, color = "CountSketch",linetype = 'CountSketch'), linewidth = 2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1,color = 'iid'), width = 0.2*x[1]) +
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2,color='Hadamard'), width = 0.2*x[1]) +
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3,color='CountSketch'), width = 0.2*x[1]) +
  scale_shape_manual(values = c("iid" = 1, "Hadamard" = 2, "CountSketch" = 3)) +
  scale_color_manual(values = c("iid" = 'green', "Hadamard" = 'red', "CountSketch" = 'blue')) +
  scale_linetype_manual(values=c("iid" = 'solid', "Hadamard" = 'longdash', "CountSketch" = 'twodash'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_gray()+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.key.height = unit(2, "line"),
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(color = "grey20", size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(color = "grey20", size = 20)
  ) +
  xlab('m') +
  ylab('coverage')+
  ggtitle('Complete sketching')

even<-seq(2,ncol(conf_pa),2);odd<-seq(1,ncol(conf_pa)-1,2)
right<-conf_pa[,even];left<-conf_pa[,odd]
accept<-(left<ls[1])&(ls[1]<right)
cov0<-apply(accept,2,sum)
cov<-matrix(cov0,ncol(conf_pa)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio<-cov/sim  #  rows represent different methods, and columns represent different sketching size m
ratio
ci<-function(x){exactci(x,sim,0.95)$conf_pa.int[1:2]} 
int=matrix(0,ncol(conf_pa)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf_pa)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    int[i,(2*j-1):(2*j)]=round(ci(cov[i,j]),digits=3)
  }
}

lower<-int[,seq(1,2*length(grid_m),2)]
upper<-int[,seq(2,2*length(grid_m),2)]


x <- grid_m
df_full <- data.frame(x = grid_m, y1 =ratio[1,], y2=ratio[2,],y3=ratio[3,], y_upper1 =upper[1,] ,
                      y_lower1 = lower[1,],y_upper2=upper[2,],y_lower2=lower[2,],y_upper3=upper[3,],y_lower3=lower[3,])

p1f_pa<-ggplot(df_full, aes(x = x)) +
  ylim(0.9, 0.98) +
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 2) +
  geom_point(aes(y = y2, color = "Hadamard",shape="Hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "Hadamard",linetype = 'Hadamard'), linewidth = 2) +
  geom_point(aes(y = y3, color = "CountSketch",shape="CountSketch"), size = 2) +
  geom_line(aes(y = y3, color = "CountSketch",linetype = 'CountSketch'), linewidth = 2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1,color = 'iid'), width = 0.2*x[1]) +
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2,color='Hadamard'), width = 0.2*x[1]) +
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3,color='CountSketch'), width = 0.2*x[1]) +
  scale_shape_manual(values = c("iid" = 1, "Hadamard" = 2, "CountSketch" = 3)) +
  scale_color_manual(values = c("iid" = 'green', "Hadamard" = 'red', "CountSketch" = 'blue')) +
  scale_linetype_manual(values=c("iid" = 'solid', "Hadamard" = 'longdash', "CountSketch" = 'twodash'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_gray()+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.key.height = unit(2, "line"),
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(color = "grey20", size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(color = "grey20", size = 20)
  ) +
  xlab('m') +
  ylab('coverage')+
  ggtitle('Partial sketching')
