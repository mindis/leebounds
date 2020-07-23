rm(list=ls())

library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))
library(geomnet)
library(ggforce)
aspect_ratio<-1.5
my_path<-"/net/holyparkesec/data/tata/leebounds/"
Figure5_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure5_selected_covs.csv"))
Figure5_dataset$group<-as.factor(Figure5_dataset$group)
ggplot(data=Figure5_dataset)+aes(x=weeks,y=bound,group=group)+
  xlab("Weeks since random assignment")+
  geom_point(aes(col=group),size=3,shape = 21,fill = "white")+
  geom_smooth(aes(col=group),se=F)+
  ylab("Expected  log wage in control status ")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c("gray","black"),
                     labels = c("lower bound", "upper bound"))+ theme(legend.position = "none")+ 
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks=c(1.2,1.6,1.8,2.0))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure5_selected_covs.png"),height=7,width=aspect_ratio*7)