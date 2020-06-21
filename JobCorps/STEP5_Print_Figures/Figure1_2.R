rm(list=ls())

library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))


my_path<-"/net/holyparkesec/data/tata/leebounds/"
aspect_ratio<-1.5


#### REPLICATE FIGURES 1 and 2
Figure1_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP3_MISC/csv/Figure1_selected_covs.csv"))
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
ggsave(paste0(my_path,"/JobCorps/Figures/Figure1_selected.png"),height=7,width=aspect_ratio*7)


Figure2_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure2_selected_covs.csv"))
test_result<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/test_result.csv"))

weeks<-1:208
weeks_001<-test_result$pvalue<=0.01
weeks_05<-test_result$pvalue<=0.05

Figure2_dataset$pvalue<-test_result$pvalue
ggplot()+
  geom_rect(data = data.frame(xmin = 91,
                              xmax = 102,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 115,
                              xmax = 117,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 121,
                              xmax = 135,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 168,
                              xmax = 173,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 177,
                              xmax = 181,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 194,
                              xmax = 208,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5)+
  geom_rect(data = data.frame(xmin = 64,
                              xmax = 68,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.3)+
  geom_rect(data = data.frame(xmin = 74,
                              xmax = 77,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.3)+
  geom_rect(data = data.frame(xmin = 102,
                              xmax = 115,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_rect(data = data.frame(xmin = 117,
                              xmax = 121,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_rect(data = data.frame(xmin = 135,
                              xmax = 156,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_rect(data = data.frame(xmin = 158,
                              xmax = 172,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_rect(data = data.frame(xmin = 173,
                              xmax = 177,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_rect(data = data.frame(xmin = 178,
                              xmax = 194,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.2)+
  geom_line(data=Figure2_dataset,
            aes(x=weeks,y=fraction),
            lwd=1)+xlab("Weeks since random assignment")+
  ylab(TeX("Fraction of subjects in covariate group $\\Delta(x)>0$"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
   scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75),expand=c(0,0),limits=c(0,0.9))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure2_selected.png"),height=7,width=aspect_ratio*7)

