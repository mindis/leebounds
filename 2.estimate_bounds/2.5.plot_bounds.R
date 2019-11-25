my_path<-"/net/holyparkesec/data/tata/Vira/"
library(feather)
library(quantreg)
library(expm)
library(purrr)
library(hdm)

selection_function_name="glm"
quantile_grid_size="0.01"
setwd(paste0(my_path,"/R_code/2.estimate_bounds"))
source(paste0(my_path,"/R_code/leebounds.R"))
#source(paste0(my_path,"/R_code/libraries.R"))
source(paste0(my_path,"/R_code/ortholeebounds.R"))
source(paste0(my_path,"/R_code/utils.R"))

estimated_bounds_CI_lee<-data.frame(test=rep(0,4))
estimated_orthobounds_CI<-data.frame(test=rep(0,4))
for (year in 1:4) {
  min_week<-(year-1)*52+1
  max_week<-year*52
  res<-read.table(paste0("Estimated_Bounds/estimated_basic_leebounds_CI_weeks_",min_week,"_",max_week,".csv"),
                                      sep=",")
  estimated_bounds_CI_lee<-cbind(estimated_bounds_CI_lee,res)
  res_ortho<-read.table(paste0("Estimated_Bounds/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"),
                        sep=",")
  res_ortho<-as.matrix(res_ortho[2:5,2:53],nrow=4,ncol=52)
  estimated_orthobounds_CI<-cbind(estimated_orthobounds_CI,res_ortho)
  
  
}
estimated_orthobounds_CI_ortho= subset(estimated_orthobounds_CI, select=-c(test))
estimated_bounds_CI_lee=subset(estimated_bounds_CI_lee,select=-c(test))

estimated_leebounds<-estimated_bounds_CI_lee[1:2,]
estimated_leebounds_CI<-estimated_bounds_CI_lee[3:4,]


estimated_orthobounds<-as.matrix(estimated_orthobounds_CI_ortho[1:2,])
estimated_orthobounds_CI<-as.matrix(estimated_orthobounds_CI_ortho[3:4,])



plotname=paste0("Figures/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),".png")
weeks=1:208                
png(plotname,width = 880,height=540)
plot(weeks, estimated_leebounds[1,],type="l",lwd=2
     ,xlab="Week since random assignment",ylab="Treatment effect",col="blue",
     ylim=c(-0.4,0.3),cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)+
  title("Lee bounds (blue) and orthogonal bounds (red) on ATE")+
  lines(weeks,estimated_orthobounds[1,],col="red",lwd=2,cex=1)+
  lines(weeks,estimated_orthobounds[2,],col="red",lwd=2,cex=1)+
  lines(weeks,estimated_leebounds[2,],col="blue",lwd=2,cex=1)+
  lines(weeks,estimated_orthobounds_CI[1,],type="l",col="red",lty=2, lwd=3)+
  lines(weeks,estimated_orthobounds_CI[2,],type="l",col="red",lty=2, lwd=3)+
  lines(weeks,estimated_leebounds_CI[1,],type="l",col="blue",lty=2, lwd=3)+
  lines(weeks,estimated_leebounds_CI[2,],type="l",col="blue",lty=2, lwd=3)+
  abline(a=0,b=0)
  legend(140,-0.2,  legend=c("Ortho bound, all covs", "Ortho bound, CI","Basic bound, no covs", "Basic bound, CI"),
       col=c("red", "red","blue","blue"), lty=1:2, cex=1.00,bty="n", lwd=3)
dev.off()






selection_function_name="rlassologit"
quantile_grid_size="0.01"

estimated_bounds_CI_lee<-data.frame(test=rep(0,4))
estimated_orthobounds_CI<-data.frame(test=rep(0,4))
estimated_leebounds_nm<-data.frame(test=rep(0,4))
for (year in 1:4) {
  min_week<-(year-1)*52+1
  max_week<-year*52
  res<-read.table(paste0("Estimated_Bounds/estimated_basic_leebounds_CI_weeks_",min_week,"_",max_week,".csv"),
                  sep=",")
  estimated_bounds_CI_lee<-cbind(estimated_bounds_CI_lee,res)
  res_ortho<-read.table(paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"),
                        sep=",")
  res_ortho<-as.matrix(res_ortho[2:5,2:53],nrow=4,ncol=52)
  estimated_orthobounds_CI<-cbind(estimated_orthobounds_CI,res_ortho)
  
  res_basic_nm<-read.table(paste0("Estimated_Bounds/estimated_leebounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"),
                           sep=",")
  res_basic_nm<-as.matrix(res_basic_nm[2:5,2:53],nrow=4,ncol=52)
  estimated_leebounds_nm<-cbind(estimated_leebounds_nm,res_basic_nm)
  
}
estimated_orthobounds_CI_ortho= subset(estimated_orthobounds_CI, select=-c(test))
estimated_bounds_CI_lee=subset(estimated_bounds_CI_lee,select=-c(test))
estimated_leebounds_nm=subset(estimated_leebounds_nm,select=-c(test))
  
estimated_leebounds<-estimated_bounds_CI_lee[1:2,]
estimated_leebounds_CI<-estimated_bounds_CI_lee[3:4,]


estimated_orthobounds<-as.matrix(estimated_orthobounds_CI_ortho[1:2,])
estimated_orthobounds_CI<-as.matrix(estimated_orthobounds_CI_ortho[3:4,])



plotname=paste0("Figures/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),".png")
weeks=1:208                
png(plotname,width = 880,height=540)
plot(weeks, estimated_leebounds_nm[1,],type="l",lwd=2
     ,xlab="Week since random assignment",ylab="Treatment effect",col="blue",
     ylim=c(-0.4,0.3),cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)+
  title("Lee bounds (blue) and orthogonal bounds (red) on ATE")+
  lines(weeks,estimated_orthobounds[1,],col="red",lwd=2,cex=1)+
  lines(weeks,estimated_orthobounds[2,],col="red",lwd=2,cex=1)+
  lines(weeks,estimated_leebounds_nm[2,],col="blue",lwd=2,cex=1)+
  lines(weeks,estimated_orthobounds_CI[1,],type="l",col="red",lty=2, lwd=3)+
  lines(weeks,estimated_orthobounds_CI[2,],type="l",col="red",lty=2, lwd=3)+
  lines(weeks,estimated_leebounds[1,],type="l",col="black",lty=2, lwd=3)+
  lines(weeks,estimated_leebounds[2,],type="l",col="black",lty=2, lwd=3)+
  abline(a=0,b=0)
legend(130,-0.2,  legend=c("Ortho bound, all covs", "Ortho bound, CI","Basic bound, no covs", "Basic bound, CI"),
       col=c("red", "red","blue","blue"), lty=1:2, cex=1.00,bty="n", lwd=3)
dev.off()



if (FALSE) {
  plot_lee_orthobounds_with_pointwise_CI(estimated_orthobounds=estimated_orthobounds,
                                         estimated_leebounds=estimated_leebounds,
                                         estimated_leebounds_CI=estimated_leebounds_CI,
                                         estimated_orthobounds_CI=estimated_orthobounds_CI,
                                         ,ylims=c(-0.2,0.2),year=floor(max_week/52))
  
  
  
}

