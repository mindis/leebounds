library(ggplot2)
library(ggforce)
library(latex2exp)
### see derivation of R, R_pointwise, beta10, beta20 in the log file 
R=0.03493337
R_pointwise=0.08823765
beta10=0.028
beta20=0.026
beta10_pointwise=0.028
beta20_pointwise=0.028

dat=data.frame(x=c(beta10,beta10_pointwise),y=c(beta20,beta20_pointwise))

df=data.frame(x=c(-0.125/2,-0.0510/2,0.02575,0.065),y=c(0.125/2,0.0510/2,-0.02575,-0.065))

ggplot(df) + 
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_x_continuous(breaks = c(-0.2,-0.1, 0,0.1,0.2),limits=c(-0.2,0.2))+
  scale_y_continuous(breaks=c(-0.2,-0.1, 0,0.1,0.2),limits=c(-0.2,0.2))+
  xlab(TeX(sprintf('Effect in week 104 ($\\beta_{104}$)')))+
  ylab(TeX(sprintf('Effect in week 208 ($\\beta_{208}$)')))+
  geom_abline(slope=-1,intercept=0,linetype=1,col="gray")+
  geom_abline(slope=1,intercept=0,linetype=1,col="gray")+
  geom_abline(slope=1,intercept=-0.128,col="black",linetype="dashed")+
  geom_abline(slope=1,intercept=-0.0515,col="black",linetype="dashed")+
  geom_abline(slope=1,intercept=0.05,col="black",linetype="dashed")+
  geom_abline(slope=1,intercept=0.125,col="black",linetype="dashed")+
  geom_point(aes(x=x,y=y),col="black")+
  annotate("text", x=-0.175, y=0.185, label=TeX(sprintf('$\\beta_{104}+\\beta_{208} = 0 $')),size=5,angle=-45)+
  annotate("text", x=-0.150, y=-0.090, label=TeX(sprintf('$-\\beta_{104}+\\beta_{208} = 0.05 $')),size=6.5,angle=45)+
  #annotate("text", x=-0.131, y=-0.000, label=TeX(sprintf('$-\\beta_{104}+\\beta_{208} = 0.13 $')),size=8,angle=45)+
  annotate("text", x=0.175, y=0.185, label=TeX(sprintf('$-\\beta_{104}+\\beta_{208} = 0 $')),size=5,angle=45)+
  annotate("text", x=-0.095, y=-0.135, label=TeX(sprintf('$-\\beta_{104}+\\beta_{208} = -0.05 $')),size=6.5,angle=45)+
  #annotate("text", x=-0.05, y=-0.17, label=TeX(sprintf('$-\\beta_{104}+\\beta_{208} = -0.13 $')),size=8,angle=45)+
  geom_circle(aes(x0 = beta10, y0 = beta20, r = R),linetype=1,col="red",lwd=1.1)+
  geom_circle(aes(x0 = beta10_pointwise, y0 = beta20_pointwise, r = R_pointwise),linetype="dotdash",col="red",lwd=1.1)
ggsave(paste0(my_path,"/JobCorps/Figures/Figure4_selected.png"),height=7,width=7)