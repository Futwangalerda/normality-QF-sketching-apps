### Figure 5

library(ggplot2)
library(tidyr)
library(patchwork)
setwd('C:/Users/Desktop/test/time')
alpha=0.1
data_time<-read.csv("C:/Users/Desktop/test/time/sum.csv")
y1_lower = data_time$srht- qnorm(1-alpha/2,sd=sqrt(data_time$srhtvar)) 
y1_upper = data_time$srht- qnorm(alpha/2,sd=sqrt(data_time$srhtvar))
y2_lower = data_time$countsketch- qnorm(1-alpha/2,sd=sqrt(data_time$countsketchvar)) 
y2_upper = data_time$countsketch- qnorm(alpha/2,sd=sqrt(data_time$countsketchvar))
y3_lower = data_time$sse- qnorm(1-alpha/2,sd=sqrt(data_time$ssevar)) 
y3_upper = data_time$sse- qnorm(alpha/2,sd=sqrt(data_time$ssevar))
y4_lower = data_time$iid- qnorm(1-alpha/2,sd=sqrt(data_time$iidvar)) 
y4_upper = data_time$iid- qnorm(alpha/2,sd=sqrt(data_time$iidvar))
y5_lower = data_time$haar- qnorm(1-alpha/2,sd=sqrt(data_time$haarvar)) 
y5_upper = data_time$haar- qnorm(alpha/2,sd=sqrt(data_time$haarvar))

data_time1<-data.frame(x = data_time$sketch_size, y1=data_time$srht, y2=data_time$countsketch, y3=data_time$sse, y4=data_time$iid, y5=data_time$haar, y1_lower=y1_lower, y1_upper=y1_upper, y2_lower=y2_lower, y2_upper=y2_upper, y3_lower=y3_lower, y3_upper=y3_upper, y4_lower=y4_lower, y4_upper=y4_upper, y5_lower=y5_lower, y5_upper=y5_upper)

p2f<-ggplot(data_time1,aes(x = x))+
  geom_point(aes(y = y1, color = "Hadamard",shape="Hadamard"), size = 2) +
  geom_line(aes(y = y1, color = "Hadamard",linetype = 'Hadamard'), linewidth = 1.2) +
  geom_point(aes(y = y2, color = "CountSketch",shape="CountSketch"), size = 2) +
  geom_line(aes(y = y2, color = "CountSketch",linetype = 'CountSketch'), linewidth = 1.2) +
  geom_point(aes(y = y3, color = "SparseSign",shape="SparseSign"), size = 2) +
  geom_line(aes(y = y3, color = "SparseSign",linetype = 'SparseSign'), linewidth = 1.2) +
  geom_point(aes(y = y4, color = "iid",shape="iid"), size = 2) +
  geom_line(aes(y = y4, color = "iid",linetype = 'iid'), linewidth = 1.2) +
  geom_point(aes(y = y5, color = "Haar",shape="Haar"), size = 2) +
  geom_line(aes(y = y5, color = "Haar",linetype = 'Haar'), linewidth = 1.2) +
  geom_errorbar(aes(ymin = y1_lower, ymax = y1_upper, color='Hadamard'), width = 50) +
  geom_errorbar(aes(ymin = y2_lower, ymax = y2_upper, color='CountSketch'), width = 50) +
  geom_errorbar(aes(ymin = y3_lower, ymax = y3_upper, color='SparseSign'), width = 50) +
  geom_errorbar(aes(ymin = y4_lower, ymax = y4_upper, color='iid'), width = 50) +
  geom_errorbar(aes(ymin = y5_lower, ymax = y5_upper, color='Haar'), width = 50) +
  scale_shape_manual(values = c("Hadamard" = 1, "CountSketch" = 2, "SparseSign" = 3, "iid" = 4, "Haar"=5)) +
  scale_color_manual(values = c("Hadamard" = 'red', "CountSketch" = 'blue', "SparseSign" = 'orange', "iid" = 'green', "Haar"='purple')) +
  scale_linetype_manual(values=c("Hadamard" = 'longdash', "CountSketch" = 'twodash', "SparseSign" = 'dotdash', "iid" = 'solid', "Haar"='dotted'))+
  guides(shape = guide_legend(title = "Method"),color = guide_legend('Method'),linetype=guide_legend('Method')) +
  labs(shape = "Merged legend",colour = "Merged legend")+
  theme_gray()+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.key.height = unit(2, "line"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 19, hjust = 0.5),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(color = "grey20", size = 20),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(color = "grey20", size = 20)
  ) +
  xlab('m') +
  ylab('Running Time in Seconds')+
  ggtitle('Sketching SVD')

p2f2 <- p2f + xlim(150,1650) + ylim(0,0.0052) + theme_gray() + theme(legend.position = 'none') +
  xlab("") + ylab("") + ggtitle("")
p2f+
  geom_rect(aes(xmin = 150, xmax = 1650, ymin = -0.01, ymax = 0.01),
            fill = "transparent", color = "black", alpha = 0, 
            linetype = "solid", linewidth = 1) +
  geom_segment(aes(x = 900, xend = 610, y = 0.01, yend = 0.29), 
               col = "black", linewidth = 0.4, linetype = "solid",
               arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  inset_element(p2f2, 0.01, 0.35, 0.6, 0.99, on_top = TRUE)

file_path<-"C:/Users/Desktop/test/time/"
ggsave("testtime_inset_1.eps", path=file_path, width = 8.5, height = 6, units = "in", dpi = 300)

