
### read the functions and also HGDP data
library(ggplot2)
library(tidyr)
library(patchwork)
library(mvtnorm)
library(phangorn)
library(PropCIs)

setwd('C:/Users/Desktop/test')
source("sketching_methods.R")
source("pivotal_inference.R")

file_path <- "C:/Users/Desktop/test/data/hgdp.txt"
file_path_save <- "C:/Users/Desktop/test/plots/1_appendix"

hgdp<-read.table(file_path,sep=',')

####  select the first 5 features to form X, and the next feature as y
X0<-as.matrix(hgdp[,1:5])
n<-nrow(X0);p<-ncol(X0);k=1;alpha=0.1
c<-c(rep(0,p-1),1)
X<-scale(X0)
n=nrow(X);p=ncol(X)
n<-2^floor(log2(n))
X<-scale(X[1:n,])/sqrt(n)
k=1;alpha=0.1
grid_m=seq(400,800,50)
lambda<-svd(X)$d
RightSingVecs<-svd(X)$v
index<-order(lambda,decreasing = TRUE)
sigma<-lambda[index]^2
RSV<-RightSingVecs[,index]
lambdaj<-sigma[k]
Projvecj<-abs(sum(c*RSV[,k]))
ci<-function(x){exactci(x,sim,0.95)$conf.int[1:2]}
ev1<-sigma[k]
ev2<-abs(sum(c*RSV[,k]))

#################################
### simulations
# sketch-and-solve for i.i.d., Hadamard, and Countsketch

sim=500
th_v_hadam<-3*(1-grid_m/n)*lambdaj^2
th_v_iid<-rep(2*lambdaj^2,length(grid_m))
th_v_countsketch<-th_v_iid

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

v1[["iid_cs_theory"]]<-th_v_iid
v1[["hadamard_theory"]]<-th_v_hadam
res_vr1<-gather(v1,type,var,-sketch_size)

data_val<-data.frame(x = v1$sketch_size, y1=v1$iid, y2=v1$hadamard, y3=v1$countsketch, y4=v1$iid_cs_theory, y5=v1$hadamard_theory)

df_var2<-data.frame(sketch_size,
                   iid=as.vector(projvec_iid)*sqrt(sketch_size),hadamard=as.vector(projvec_srht)*sqrt(sketch_size),countsketch=as.vector(projvec_countsketch)*sqrt(sketch_size))

v2<-aggregate(.~sketch_size,df_var2, FUN = stats::var)

v2[["iid_cs_theory"]]<-th_vec_iid
v2[["hadamard_theory"]]<-th_vec_hadam

data_vec<-data.frame(x = v2$sketch_size, y1=v2$iid, y2=v2$hadamard, y3=v2$countsketch, y4=v2$iid_cs_theory, y5=v2$hadamard_theory)


p1f<-ggplot(data_val,aes(x = x))+
 geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
 geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
 geom_point(aes(y = y2, color = "hadamard",shape="hadamard"), size = 2) +
 geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
 geom_point(aes(y = y3, color = "countsketch",shape="countsketch"), size = 2) +
 geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
 geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 2) +
 geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.2) +
 geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 2) +
 geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.2) +
 geom_errorbar(aes(ymin= y1 -qnorm(1-alpha/2,sd=sqrt(2*y1^2/(x-1))), ymax=y1-qnorm(alpha/2,sd=sqrt(2*y1^2/(x-1))), color="iid"), width=10)+
 geom_errorbar(aes(ymin= y2 -qnorm(1-alpha/2,sd=sqrt(2*y2^2/(x-1))), ymax=y2-qnorm(alpha/2,sd=sqrt(2*y2^2/(x-1))), color="hadamard"), width=10)+
 geom_errorbar(aes(ymin= y3 -qnorm(1-alpha/2,sd=sqrt(2*y3^2/(x-1))), ymax=y3-qnorm(alpha/2,sd=sqrt(2*y3^2/(x-1))), color="countsketch"), width=10)+
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
  ylab(NULL)
  #ylab('Coverage')+
  #ggtitle('Eigenvalue')
ggsave("k1_HGDP_varval_new.eps", path=file_path_save, width = 4, height = 6, units = "in", dpi = 300)


p2f<-ggplot(data_vec,aes(x = x))+
 geom_point(aes(y = y1, color = "iid",shape='iid'), size = 2) +
 geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
 geom_point(aes(y = y2, color = "hadamard",shape="hadamard"), size = 2) +
 geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
 geom_point(aes(y = y3, color = "countsketch",shape="countsketch"), size = 2) +
 geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
 geom_point(aes(y = y4, color = "iid_cs_theory",shape='iid_cs_theory'), size = 2) +
 geom_line(aes(y = y4, color = "iid_cs_theory",linetype = 'iid_cs_theory'), linewidth = 1.2) +
 geom_point(aes(y = y5, color = "hadamard_theory",shape="hadamard_theory"), size = 2) +
 geom_line(aes(y = y5, color = "hadamard_theory",linetype = 'hadamard_theory'), linewidth = 1.2) +
 geom_errorbar(aes(ymin= y1 -qnorm(1-alpha/2,sd=sqrt(2*y1^2/(x-1))), ymax=y1-qnorm(alpha/2,sd=sqrt(2*y1^2/(x-1))), color="iid"), width=10)+
 geom_errorbar(aes(ymin= y2 -qnorm(1-alpha/2,sd=sqrt(2*y2^2/(x-1))), ymax=y2-qnorm(alpha/2,sd=sqrt(2*y2^2/(x-1))), color="hadamard"), width=10)+
 geom_errorbar(aes(ymin= y3 -qnorm(1-alpha/2,sd=sqrt(2*y3^2/(x-1))), ymax=y3-qnorm(alpha/2,sd=sqrt(2*y3^2/(x-1))), color="countsketch"), width=10)+
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
  ylab(NULL)
  #ylab('Coverage')+
  #ggtitle('Eigenvector')
ggsave("k1_HGDP_varvec_new.eps", path=file_path_save, width = 4, height = 6, units = "in", dpi = 300)

#################################################
### comparison of the coverage

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
data_val_full<-data.frame(x=grid_m, y1 =ratio1[1,], y2=ratio1[2,],y3=ratio1[3,], y_upper1 =upper1[1,],
                        y_lower1 = lower1[1,],y_upper2=upper1[2,],y_lower2=lower1[2,],y_upper3=upper1[3,],y_lower3=lower1[3,])

intvec=matrix(0,ncol(conf2)/(2*length(grid_m)),2*length(grid_m))
for(i in 1:(ncol(conf2)/(2*length(grid_m)))){
  for(j in 1:length(grid_m)){
    intvec[i,(2*j-1):(2*j)]=round(ci(covar2[i,j]),digits=3)
  }
}

lower2<-intvec[,seq(1,2*length(grid_m),2)]
upper2<-intvec[,seq(2,2*length(grid_m),2)]

data_vec_full<-data.frame(x=grid_m, y1 =ratio2[1,], y2=ratio2[2,],y3=ratio2[3,], y_upper1 =upper2[1,],
                        y_lower1 = lower2[1,],y_upper2=upper2[2,],y_lower2=lower2[2,],y_upper3=upper2[3,],y_lower3=lower2[3,])

# write.csv(df_val_full,"C:/Users/Desktop/test/plots/csv/k1_HGDP_Coverage_val.csv")
# write.csv(df_vec_full,"C:/Users/Desktop/test/plots/csv/k1_HGDP_Coverage_vec.csv")

p1f<-ggplot(data_val_full,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape= 'iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "hadamard",shape= "hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "countsketch",shape= "countsketch"), size = 2) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1, color = 'iid'), width=0.05*data_val_full$x[1])+
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2, color = "hadamard"), width=0.05*data_val_full$x[1])+
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3, color = "countsketch"), width=0.05*data_val_full$x[1])+
  scale_shape_manual(values = c("iid" = 1, "hadamard" = 2, "countsketch" = 3)) +
  scale_color_manual(values = c("iid" = 'green', "hadamard" = 'red', "countsketch" = 'blue')) +
  scale_linetype_manual(values=c("iid" = 'solid', "hadamard" = 'longdash', "countsketch" = 'twodash'))+
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
  ylab(NULL)
  #ylab('Coverage')+
  #ggtitle('Eigenvalue')
ggsave("k1_HGDP_Coverage_val_new.eps", path="C:/Users/Desktop/test/plots/1_appendix", width = 4, height = 6, units = "in", dpi = 300)

p2f<-ggplot(data_vec_full,aes(x = x))+
  geom_point(aes(y = y1, color = "iid",shape= 'iid'), size = 2) +
  geom_line(aes(y = y1, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "hadamard",shape= "hadamard"), size = 2) +
  geom_line(aes(y = y2, color = "hadamard",linetype = 'hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "countsketch",shape= "countsketch"), size = 2) +
  geom_line(aes(y = y3, color = "countsketch",linetype = 'countsketch'), linewidth = 1.2) +
  geom_errorbar(aes(ymin = y_lower1, ymax = y_upper1, color = 'iid'), width=0.05*data_vec_full$x[1])+
  geom_errorbar(aes(ymin = y_lower2, ymax = y_upper2, color = "hadamard"), width=0.05*data_vec_full$x[1])+
  geom_errorbar(aes(ymin = y_lower3, ymax = y_upper3, color = "countsketch"), width=0.05*data_vec_full$x[1])+
  scale_shape_manual(values = c("iid" = 1, "hadamard" = 2, "countsketch" = 3)) +
  scale_color_manual(values = c("iid" = 'green', "hadamard" = 'red', "countsketch" = 'blue')) +
  scale_linetype_manual(values=c("iid" = 'solid', "hadamard" = 'longdash', "countsketch" = 'twodash'))+
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
  ylab(NULL)
  #ylab('Coverage')+
  #ggtitle('Eigenvector')
ggsave("k1_HGDP_Coverage_vec_new.eps", path="C:/Users/Desktop/test/plots/1_appendix", width = 4, height = 6, units = "in", dpi = 300)

