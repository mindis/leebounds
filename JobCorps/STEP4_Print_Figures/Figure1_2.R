rm(list=ls())

library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))
library(geomnet)
library(ggforce)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
aspect_ratio<-1.5


#### REPLICATE FIGURES 1 and 2
Figure1_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP3_MISC/csv/Figure1.csv"))
Figure1_dataset$group<-as.factor(Figure1_dataset$group)


ggplot(data=Figure1_dataset)+aes(x=weeks,y=delta,group=group)+
  geom_point(aes(shape=group),size=2)+xlab("Weeks since random assignment")+
  ylab(TeX(sprintf('$\\Delta=Prob (S=1|D=1) - Prob(S=1|D=0)$')))+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0)+
  scale_shape_manual(values=c(19,1),
                     guide = guide_legend(reverse = TRUE))+ theme(legend.position = "none")+ 
  annotate("text", x=120, y=0.1, label=TeX(sprintf('$x:  \\Delta(x)>0$')),size=6)+ 
  annotate("text", x=120, y=-0.1, label=TeX(sprintf('$x:  \\Delta(x)<0$')),size=6)+
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks=c(-0.1,0.,0.1))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure1.png"),height=7,width=aspect_ratio*7)


Figure2_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP3_MISC/csv/Figure2.csv"))
ggplot(data=Figure2_dataset)+aes(x=weeks,y=fraction)+
  geom_line(lwd=1)+xlab("Weeks since random assignment")+
  ylab(TeX("Fraction of subjects in covariate group $\\Delta(x)>0$"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75),expand=c(0,0),limits=c(0,0.9))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure2.png"),height=7,width=aspect_ratio*7)

