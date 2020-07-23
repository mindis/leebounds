rm(list=ls())

library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))


my_path<-"/net/holyparkesec/data/tata/leebounds/"
aspect_ratio<-1.5


#### REPLICATE FIGURES 1 and 2
Figure1_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure1_selected_covs.csv"))
Figure1_dataset$group<-as.factor(Figure1_dataset$group)
Figure2_dataset<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure2_selected_covs.csv"))
test_result<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/test_result.csv"))


#Figure1_dataset$fraction<-sapply(c(Figure2_dataset$fraction,1-Figure2_dataset$fraction),min,0.5)
Figure1_dataset$fraction<-c(Figure2_dataset$fraction,1-Figure2_dataset$fraction)
                                 
ggplot(data=Figure1_dataset)+aes(x=weeks,y=delta,group=group,size=fraction/2)+
  geom_point(aes(color=group))+xlab("Weeks since random assignment")+
  ylab("Treatment-control difference in employment rate")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0)+
  scale_color_manual(values=c("gray","black"))+ theme(legend.position = "none")+ 
  annotate("text", x=120, y=0.1, label='Predicted Positive Effect',size=6)+ 
  annotate("text", x=120, y=-0.1, label='Predicted Negative Effect',size=6)+
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks=c(-0.1,0.,0.1))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure1_selected.png"),height=7,width=aspect_ratio*7)



sum(test_result$pvalue<=0.01)
sum(test_result$pvalue>0.01 & test_result$pvalue<=0.05)
weeks<-1:208
weeks_001<-test_result$pvalue<=0.01
weeks_05<-test_result$pvalue<=0.05

Figure2_dataset$pvalue<-test_result$pvalue
Figure2_dataset$group<-"pvalue > 0.05"
Figure2_dataset$group[Figure2_dataset$pvalue <=0.01]<-"pvalue < 0.01"
Figure2_dataset$group[Figure2_dataset$pvalue >0.01 & Figure2_dataset$pvalue<=0.05 ]<-"pvalue in [0.01, 0.05]"

ggplot(data=Figure2_dataset)+
  geom_point( aes(x=weeks,y=pvalue,color=group))+
  geom_point( aes(x=weeks,y=pvalue),color="white")+
  scale_color_manual(values=c("gray100","gray90","gray75"),name="Monotonicity test p-value",
                     labels=c("> 0.05", "[0.01, 0.05]", "< 0.01" ))+
  geom_rect(data = data.frame(xmin = 91,
                              xmax = 102,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 115,
                              xmax = 117,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 121,
                              xmax = 135,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 168,
                              xmax = 173,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 177,
                              xmax = 181,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 194,
                              xmax = 208,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 64,
                              xmax = 68,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 74,
                              xmax = 77,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 102,
                              xmax = 115,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 117,
                              xmax = 121,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 135,
                              xmax = 156,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 158,
                              xmax = 172,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 173,
                              xmax = 177,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 178,
                              xmax = 194,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_line(data=Figure2_dataset,
            aes(x=weeks,y=fraction),
            lwd=1)+xlab("Weeks since random assignment")+
  ylab(TeX("Fraction of applicants with positive employment effect"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+theme(legend.position="bottom",legend.text=element_text(size=15),legend.title = element_text(size=15))
   scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75),expand=c(0,0),limits=c(0,0.9))

ggsave(paste0(my_path,"/JobCorps/Figures/Figure2_selected.png"),height=7,width=aspect_ratio*7)

