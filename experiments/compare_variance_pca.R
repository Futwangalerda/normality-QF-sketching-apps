### Figure 4, 9
### Empirical Asymptotic Variance Comparison

library(mvtnorm)
library(phangorn)
library(ggplot2)
library(tidyr)
setwd('C:/Users/Desktop/test')
source("sketching_methods.R")
file_path<-"C:/Users/Desktop/test/plots"

##Estimating the largest singular value
n=2^11;p=15;alpha=0.1;k=1
grid_m=seq(200,1600,200) # 200-1600, 200
c<-c(1,rep(0,p-1))
set.seed(15)
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
#   A<-matrix(0,n,n)
#   for(i in 1:n){
#     for(j in 1:n){
#       A[i,j]=t^(abs(i-j))
#     }
#   }
#   return(A)
# }
# sigX<-2*genesig(0.5,p)
# W<-rmvt(n, sigX, df = 8)
# U<-svd(W)$u
# V<-svd(matrix(rnorm(p*p),p))$v
# D<-diag(seq(0.1,1,by=(1-0.1)/(p-1)))
# X<-U%*%D%*%V

### Case 3
# X1<-matrix(rnorm(n/2*p),n/2,p)
# X2<-matrix(rnorm(n/2*p,5,1),n/2,p)
# X<-rbind(X1,X2)
# X<-scale(X)/sqrt(n)

lambda<-svd(X)$d
RightSingVecs<-svd(X)$v
index<-order(lambda,decreasing = TRUE)
sigma<-lambda[index]^2
RSV<-RightSingVecs[,index]
lambdaj<-sigma[k]
Projvecj<-abs(sum(c*RSV[,k]))

####################
##theoretical asymptotic variances for Hadamard sketching
th_val_hadam<-3*(1-grid_m/n)*lambdaj^2 

####################
###theoretical asymptotic variances for i.i.d. sketching estimators and Countsketch
th_val_iid<-rep(2*lambdaj^2,length(grid_m))
th_val_countsketch<-th_val_iid

####################
###theoretical asymptotic variances for Haar sketching estimators
th_val_haar<-2*(1-grid_m/n)*lambdaj^2 


###variance of the projection of the right singular vector onto one specific direction
ProjVar<-0
for(i in 1:p){
  if(i!=k){
    ProjVar<-ProjVar+sigma[i]*sigma[k]/(sigma[i]-sigma[k])^2*sum(c*RSV[,i])^2
  }
}

th_vec_hadam<-(1-grid_m/n)*ProjVar
th_vec_haar<-th_vec_hadam
th_vec_iid<-rep(ProjVar,length(grid_m))
th_vec_countsketch<-th_vec_iid

#####################
### simulations
# sketch-and-solve for i.i.d., Countsketch, and Hadamard
sim=500

r_iid<-matrix(0,sim,length(grid_m))
projvec_iid<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_iid(k,grid_m[j],X,c,type=1)
    r_iid[i,j]<-Temp$lambdai
    projvec_iid[i,j]<-Temp$squareprojveci
  }
}


r_srht<-matrix(0,sim,length(grid_m))
projvec_srht<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_SRHT(k,grid_m[j],X,c)
    r_srht[i,j]<-Temp$lambdai
    projvec_srht[i,j]<-Temp$squareprojveci
  }
}

r_countsketch<-matrix(0,sim,length(grid_m))
projvec_countsketch<-matrix(0,sim,length(grid_m))
for(i in 1:sim){
  for(j in 1:length(grid_m)){
    Temp<-Estising_CountSketch(k,grid_m[j],X,c)
    r_countsketch[i,j]<-Temp$lambdai
    projvec_countsketch[i,j]<-Temp$squareprojveci
  }
}

#################################################
### display the variances
sketch_size<-as.vector(t(matrix(rep(grid_m,sim),length(grid_m))))
df_var1<-data.frame(sketch_size,
                    iid=as.vector(r_iid)*sqrt(sketch_size),hadamard=as.vector(r_srht)*sqrt(sketch_size),countsketch=as.vector(r_countsketch)*sqrt(sketch_size))

v1<-aggregate(.~sketch_size,df_var1, FUN = stats::var)

v1[["iid_cs_theory"]]<-th_val_iid
v1[["hadamard_theory"]]<-th_val_hadam

data_val<-data.frame(x = v1$sketch_size/n, y1=v1$iid, y2=v1$hadamard, y3=v1$countsketch, y4=v1$iid_cs_theory, y5=v1$hadamard_theory, y6=v1)

df_var2<-data.frame(sketch_size,
                    iid=as.vector(projvec_iid)*sqrt(sketch_size),hadamard=as.vector(projvec_srht)*sqrt(sketch_size),countsketch=as.vector(projvec_countsketch)*sqrt(sketch_size))

v2<-aggregate(.~sketch_size,df_var2, FUN = stats::var)

v2[["iid_cs_theory"]]<-th_vec_iid
v2[["hadamard_theory"]]<-th_vec_hadam

data_vec<-data.frame(x = v2$sketch_size/n, y1=v2$iid, y2=v2$hadamard, y3=v2$countsketch, y4=v2$iid_cs_theory, y5=v2$hadamard_theory)

# write.csv(data_val,"C:/Users/Desktop/test/testplots/csv/n2048p15k1_Case2_new_val.csv")
# write.csv(data_vec,"C:/Users/Desktop/test/testplots/csv/n2048p15k1_Case2_new_vec.csv")

p1f<-ggplot(data_val,aes(x = x))+
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
  geom_errorbar(aes(ymin= y1 -qnorm(1-alpha/2,sd=sqrt(2*y1^2/(x*n-1))), ymax=y1-qnorm(alpha/2,sd=sqrt(2*y1^2/(x*n-1))), color="iid"), width=0.01)+
  geom_errorbar(aes(ymin= y2 -qnorm(1-alpha/2,sd=sqrt(2*y2^2/(x*n-1))), ymax=y2-qnorm(alpha/2,sd=sqrt(2*y2^2/(x*n-1))), color="hadamard"), width=0.01)+
  geom_errorbar(aes(ymin= y3 -qnorm(1-alpha/2,sd=sqrt(2*y3^2/(x*n-1))), ymax=y3-qnorm(alpha/2,sd=sqrt(2*y3^2/(x*n-1))), color="countsketch"), width=0.01)+
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
  ylab(NULL)+
#ylab('Coverage')+
#ylab('Eigenvalue')
  ggtitle('Eigenvalue')
ggsave("n2048p15k1_Case2_varval_new.eps", path=file_path, width = 4, height = 6, units = "in", dpi = 300)

p2f<-ggplot(data_vec,aes(x = x))+
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
  geom_errorbar(aes(ymin= y1 -qnorm(1-alpha/2,sd=sqrt(2*y1^2/(x*n-1))), ymax=y1-qnorm(alpha/2,sd=sqrt(2*y1^2/(x*n-1))), color="iid"), width=0.01)+
  geom_errorbar(aes(ymin= y2 -qnorm(1-alpha/2,sd=sqrt(2*y2^2/(x*n-1))), ymax=y2-qnorm(alpha/2,sd=sqrt(2*y2^2/(x*n-1))), color="hadamard"), width=0.01)+
  geom_errorbar(aes(ymin= y3 -qnorm(1-alpha/2,sd=sqrt(2*y3^2/(x*n-1))), ymax=y3-qnorm(alpha/2,sd=sqrt(2*y3^2/(x*n-1))), color="countsketch"), width=0.01)+
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
  ylab(NULL)+
  ggtitle('Eigenvector')
ggsave("n2048p15k1_Case2_varvec_new.eps", path=file_path, width = 4, height = 6, units = "in", dpi = 300)
