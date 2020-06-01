rm(list=ls())

my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))

### Figure 1 ###
Lee_data_summary_earn_con<-apply(Lee_data[Lee_data$TREATMNT.y==0,grep("EARNH",colnames(Lee_data),value=TRUE)],2,mean)
Lee_data_summary_earn_trt<-apply(Lee_data[Lee_data$TREATMNT.y==1,grep("EARNH",colnames(Lee_data),value=TRUE)],2,mean)

png(paste0(my_path,"/JobCorps/Figures/Lee_2009_Figure1.png"))
p<-plot(x=1:208,y=Lee_data_summary_earn_con,xlim=c(0,208),ylim=c(0,250),xlab="Week since random assignment",ylab="Weekly earnings (dollars)")+
  points(x=1:208,y=Lee_data_summary_earn_trt,xlim=c(0,208),ylim=c(0,250),pch=3)
dev.off()

Lee_data_summary_emp_con<-apply(Lee_data[Lee_data$TREATMNT.y==0,grep("HWH",colnames(Lee_data),value=TRUE)],2,function (x) mean(x>0))
Lee_data_summary_emp_trt<-apply(Lee_data[Lee_data$TREATMNT.y==1,grep("HWH",colnames(Lee_data),value=TRUE)],2,function (x) mean(x>0))


png(paste0(my_path,"/JobCorps/Figures/Lee_2009_Figure2.png"))
p<-plot(x=1:208,y=Lee_data_summary_emp_con,xlim=c(0,208),ylim=c(0,0.7),xlab="Week since random assignment",ylab="Weekly earnings (dollars)")+
  points(x=1:208,y=Lee_data_summary_emp_trt,xlim=c(0,208),ylim=c(0,0.7),pch=3)
dev.off()