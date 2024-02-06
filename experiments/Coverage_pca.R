### Figure 3, 7
### Comparison of the Coverage of Nominal 95% Intervals in Different Sketching

setwd('C:/Users/Desktop/test')
source("sketching_methods.R")
source("pivotal_inference.R")

library(mvtnorm)
library(phangorn)
library(PropCIs)
library(ggplot2)
library(tidyr)

n=2^11;p=15;k=2
sim=500
c<-c(1,rep(0,p-1))
grid_m=seq(200,1600,200)
set.seed(15)
ci<-function(x){exactci(x,sim,0.95)$conf.int[1:2]}
### Case 1
d<-1/(1:p)
D<-diag(d)
O1<-matrix(rnorm(n*p),n,p)
W<-svd(O1)$u
O2<-matrix(rnorm(p*p),p,p)
U<-qr.Q(qr(O2))
X<-W %*% D %*% t(U)

### Case 2
# genesig<-function(t,n){
# A<-matrix(0,n,n)
# for(i in 1:n){
#   for(j in 1:n){
#     A[i,j]=t^(abs(i-j))
#   }
# }
# return(A)
# }
# sigX<-2*genesig(0.5,p)
# W<-rmvt(n, sigX, df = 2)
# U<-svd(W)$u
# V<-svd(matrix(rnorm(p*p),p))$v
# D<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
# X<-U%*%D%*%V

### Case 3
# X1<-matrix(rnorm(n/2*p),n/2,p)
# X2<-matrix(rnorm(n/2*p,5,1),n/2,p)
# X<-rbind(X1,X2)

lambda<-svd(X)$d
RightSingVecs<-svd(X)$v
index<-order(lambda,decreasing = TRUE)
sigma<-lambda[index]^2
RSV<-RightSingVecs[,index]
ev1<-sigma[k]
ev2<-sqrt(sum(c*RSV[,k])^2)

pivo_conf_iid_val<-matrix(0,sim,2*length(grid_m))
pivo_conf_iid_projvec<-matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_iid(k,grid_m[j],X,c,type=1)
    SX<-Temp$matrix
    pivo_conf_iid_val[i,(2*j-1):(2*j)]<-pivo_iid_val(k,n,SX,alpha=0.05)$conf
    pivo_conf_iid_projvec[i,(2*j-1):(2*j)]<-pivo_iid_projvec(k,n,SX,alpha=0.05)$conf
  }
}

#pivo_conf_haar_val<-matrix(0,sim,2*length(grid_m))
#pivo_conf_haar_projvec<-matrix(0,sim,2*length(grid_m))
#for(i in 1:sim){
#  for(j in 1:length(grid_m)){
#    Temp<-Estising_Haar(k,grid_m[j],X,c,type=1)
#    SX<-Temp$matrix
#    pivo_conf_haar_val[i,(2*j-1):(2*j)]<-pivo_haar_val(k,n,SX,alpha=0.05)$conf
#    pivo_conf_haar_projvec[i,(2*j-1):(2*j)]<-pivo_haar_projvec(k,n,SX,alpha=0.05)$conf
#  }
#}

X<-padding(X)
n0<-nrow(X)

pivo_conf_srht_val<-matrix(0,sim,2*length(grid_m))
pivo_conf_srht_projvec<-matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_SRHT(k,grid_m[j],X,c)
    SX<-Temp$matrix
    pivo_conf_srht_val[i,(2*j-1):(2*j)]<-pivo_srht_val(k,n0,SX,alpha=0.05)$conf
    pivo_conf_srht_projvec[i,(2*j-1):(2*j)]<-pivo_srht_projvec(k,n0,SX,alpha=0.05)$conf
  }
}

pivo_conf_CountSketch_val<-matrix(0,sim,2*length(grid_m))
pivo_conf_CountSketch_projvec<-matrix(0,sim,2*length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_CountSketch(k,grid_m[j],X,c)
    SX<-Temp$matrix
    pivo_conf_CountSketch_val[i,(2*j-1):(2*j)]<-pivo_iid_val(k,n,SX,alpha=0.05)$conf
    pivo_conf_CountSketch_projvec[i,(2*j-1):(2*j)]<-pivo_iid_projvec(k,n,SX,alpha=0.05)$conf
  }
}


conf1<-cbind(pivo_conf_iid_val,pivo_conf_srht_val,pivo_conf_CountSketch_val)

even1<-seq(2,ncol(conf1),2);odd1<-seq(1,ncol(conf1)-1,2)
right1<-conf1[,even1];left1<-conf1[,odd1]
accept1<-(left1<ev1)&(ev1<right1)
cov1<-apply(accept1,2,sum)

covar1<-matrix(cov1,ncol(conf1)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio1<-covar1/sim

conf2<-cbind(pivo_conf_iid_projvec,pivo_conf_srht_projvec,pivo_conf_CountSketch_projvec)

even2<-seq(2,ncol(conf2),2);odd2<-seq(1,ncol(conf2)-1,2)
right2<-conf2[,even2];left2<-conf2[,odd2]
accept2<-(left2<ev2)&(ev2<right2)
cov2<-apply(accept2,2,sum)

covar2<-matrix(cov2,ncol(conf2)/(2*length(grid_m)),length(grid_m),byrow=TRUE)
ratio2<-covar2/sim
 
intval=matrix(0,ncol(conf1)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf1)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    intval[i,(2*j-1):(2*j)]=round(ci(covar1[i,j]),digits=3)
  }
}

lower1<-intval[,seq(1,2*length(grid_m),2)]
upper1<-intval[,seq(2,2*length(grid_m),2)]

x<-grid_m
df_val_full<-data.frame(x=grid_m, y1 =ratio1[1,], y2=ratio1[2,],y3=ratio1[3,], y_upper1 =upper1[1,],
                      y_lower1 = lower1[1,],y_upper2=upper1[2,],y_lower2=lower1[2,],y_upper3=upper1[3,],y_lower3=lower1[3,])

intvec=matrix(0,ncol(conf2)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf2)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    intvec[i,(2*j-1):(2*j)]=round(ci(covar2[i,j]),digits=3)
  }
}

lower2<-intvec[,seq(1,2*length(grid_m),2)]
upper2<-intvec[,seq(2,2*length(grid_m),2)]

df_vec_full<-data.frame(x=grid_m, y1 =ratio2[1,], y2=ratio2[2,],y3=ratio2[3,], y_upper1 =upper2[1,],
                        y_lower1 = lower2[1,],y_upper2=upper2[2,],y_lower2=lower2[2,],y_upper3=upper2[3,],y_lower3=lower2[3,])

write.csv(df_val_full,"C:/Users/Desktop/test/plots/csv/n2048p15k1_Case2_Coverage_val.csv")
write.csv(df_vec_full,"C:/Users/Desktop/test/plots/csv/n2048p15k1_Case2_Coverage_vec.csv")

p1f<-ggplot(df_val_full,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "Hadamard",shape="Hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "Hadamard",linetype = 'Hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "CountSketch",shape="CountSketch"), size = 2) +
  geom_line(aes(y = y3, color = "CountSketch",linetype = 'CountSketch'), linewidth = 1.2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1,color = 'iid'), width = 50) +
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2,color='Hadamard'), width = 50) +
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3,color='CountSketch'), width = 50) +
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
  ggtitle('Eigenvalue')

p2f<-ggplot(data_vec,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "Hadamard",shape="Hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "Hadamard",linetype = 'Hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "CountSketch",shape="CountSketch"), size = 2) +
  geom_line(aes(y = y3, color = "CountSketch",linetype = 'CountSketch'), linewidth = 1.2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1,color = 'iid'), width = 50) +
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2,color='Hadamard'), width = 50) +
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3,color='CountSketch'), width = 50) +
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
  ggtitle('Eigenvector')

