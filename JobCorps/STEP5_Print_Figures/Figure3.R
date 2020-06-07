
rm(list=ls())

library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))
library(geomnet)
library(ggforce)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
aspect_ratio<-1.5

Figure3_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure3_selected_covs.csv"))
#Figure3_dataset<-Figure3_dataset[1:110,]
ggplot(data=Figure3_dataset)+aes(x=weeks)+
  xlab("Weeks since random assignment")+
  geom_line(aes(y = fraction/5,colour="fraction"),size=0.25,key_glyph = draw_key_rect)+
  geom_point(aes(y=lb_positive,colour="lb_positive"),size=3,shape = 21,fill = "white")+
  geom_point(aes(y=ub_positive,colour="ub_positive"),size=3,shape = 21,fill = "white")+
  geom_smooth(aes(y=lb_positive,colour="lb_positive"),se=F)+
  geom_smooth(aes(y=ub_positive,colour="ub_positive"),se=F)+
  ylab(TeX('Bounds on average effect'))+
  theme_bw()+
  scale_y_continuous(sec.axis =sec_axis(~.*3, name = "Fraction of covariate groups"),
                     limits =c(0,0.33),expand=c(0,0))+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c("black","blue","red")) + theme(legend.position = "none")+ 
  theme(legend.title =element_blank(),
        legend.direction = "vertical",
        legend.text=element_text(size=15),
        legend.text.align = 0)+
  scale_x_continuous(breaks = c(40,80,120,160,200),limits=c(0,210))

ggsave(paste0(my_path,"/JobCorps/Figures/Figure3_selected_covs.png"),height=7,width=aspect_ratio*7)