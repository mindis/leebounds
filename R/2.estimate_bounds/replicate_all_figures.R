### Figure Trimming_unweighted.png and Basic_Lee_bounds_only.png are replicated in 2.1.trimming_prop_and_basic_lee_bounds.R
my_path<-"/net/holyparkesec/data/tata/leebounds/"
### Figure 2. Basic Lee bounds and 95 % CI
weeks<-1:208
# weeks from Lee (2009)
selected_weeks<-c(45,90,135,180)

basic_bounds<-read.csv(paste0(my_path,"/R/2.estimate_bounds/Estimated_Bounds/estimated_basic_leebounds_CI_weeks_1_208.csv"))
estimated_leebounds<-basic_bounds[1:2,]
estimated_leebounds_CI<-basic_bounds[3:4,]
setwd(paste0(my_path,"/R/2.estimate_bounds/"))
png("Figures/Basic_Lee_bounds.png")
plot(weeks,estimated_leebounds[2,],xlab="Weeks since random assignment",ylab="Lower and upper bound",col="blue",
     type="l",lwd=3,cex.lab=1.5, cex.axis=1.5,cex=1,ylim=c(-0.3,0.3))+
  lines(weeks,estimated_leebounds[1,],col="blue",lwd=3)+
  lines(weeks,estimated_leebounds_CI[1,],col="blue",lty=2, lwd=3)+
  lines(weeks,estimated_leebounds_CI[2,],col="blue",lty=2, lwd=3)+
  title("Lee (2009) bounds")+
  points(x=selected_weeks,y=estimated_leebounds[1,selected_weeks],col="black",lwd=5)+
  points(x=selected_weeks,y=estimated_leebounds[2,selected_weeks],col="black",lwd=5)+
  abline(a=0,b=0,col="black",lwd=2)
dev.off()

### Figure 3. Estimated proportion of people in the treatment/control group
p.0.hat.nonmonotone<-read.csv(paste0(my_path,"/R/2.estimate_bounds/First_Stage_Predicted_Values/Predicted_Selection/p.0.hat.nonmonotone.glm_weeks_1_208.csv"))
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
estimated.p.0.hat<-p.0.hat.nonmonotone[,2:209]
fraction_treat_helps<-apply(estimated.p.0.hat<1,2,mean)
png("Figures/Prop_treat_helps.png")

plot(weeks,fraction_treat_helps,xlab="Weeks since random assignment",ylab="Fraction",col="blue",
     type="l",lwd=3,cex.lab=1.5, cex.axis=1.5,cex=1,ylim=c(0,1))+
    title("Fraction of subjects whose employment JobCorps helps")
dev.off()
employment_rate_X1_treat<-rep(0,length(weeks))
employment_rate_X1_control<-rep(0,length(weeks))
employment_rate_X0_treat<-rep(0,length(weeks))
employment_rate_X0_control<-rep(0,length(weeks))


for (i in 1:length(weeks)) {
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  logwage_week[logwage_week>0]<-logwage_week[logwage_week>0]+rnorm(sum(logwage_week>0),0.001)
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  employment_rate_X0_treat[i]<-mean(leedata_week$selection[leedata_week$treat==1 & estimated.p.0.hat[,i]<=1 ],na.rm=TRUE)
  employment_rate_X0_control[i]<-mean(leedata_week$selection[leedata_week$treat==0 & estimated.p.0.hat[,i]<=1 ],na.rm=TRUE)
  
  employment_rate_X1_treat[i]<-mean(leedata_week$selection[leedata_week$treat==1 & estimated.p.0.hat[,i]>1 ],na.rm=TRUE)
  employment_rate_X1_control[i]<-mean(leedata_week$selection[leedata_week$treat==0 & estimated.p.0.hat[,i]>1 ],na.rm=TRUE)
  
  ## deciding which covariates to put in
}

png("Figures/Treat_control_diff_by_group.png")
par(mfrow=c(1,2))
plot(weeks,sapply(employment_rate_X0_treat-employment_rate_X0_control,max,0),xlab="Weeks",ylab="Treatment-control difference",col="blue",
    lwd=3,cex.lab=1.5, cex.axis=1.5,cex=1,ylim=c(-0.15,0.15))+
     title("Treatment helps")+
    abline(a=0,b=0)

plot(weeks,sapply(employment_rate_X1_treat-employment_rate_X1_control,min,0),xlab="Weeks",ylab="Treatment-control difference",col="blue",
              lwd=3,cex.lab=1.5, cex.axis=1.5,cex=1,ylim=c(-0.15,0.15))+
  title("Treatment hurts")+
  abline(a=0,b=0)
dev.off()

## Monotone bounds based on glm: wide and tight
selection_function_name<-"glm"
quantile_grid_size<-0.01
min_week<-1
max_week<-208
estimated_orthobounds_CI<-read.csv(paste0("Estimated_Bounds/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))
estimated_leebounds_CI<-read.csv(paste0("Estimated_Bounds/estimated_leebounds_monotone_weeks_",min_week,"_",max_week,".csv"))

estimated_orthobounds_CI<-estimated_orthobounds_CI[,2:209]
ylims=c(-0.4,0.3)
plotname<-paste0("Figures/monotone_",selection_function_name,".png")
plot_bounds(plotname=plotname,estimated_leebounds_CI=estimated_leebounds_CI,
            estimated_orthobounds_CI=estimated_orthobounds_CI)


estimated_orthobounds_CI<-read.csv(paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))
estimated_leebounds_CI<-read.csv(paste0("Estimated_Bounds/estimated_leebounds_nonmonotone_",selection_function_name,"_weeks_",min_week,"_",max_week,".csv"))


estimated_orthobounds_CI<-estimated_orthobounds_CI[,2:209]
estimated_leebounds_CI<-estimated_leebounds_CI[,2:209]


ylims=c(-0.4,0.3)
plotname<-paste0("Figures/nonmonotone_",selection_function_name,".png")
plot_bounds(plotname=plotname,estimated_leebounds_CI=estimated_leebounds_CI,
            estimated_orthobounds_CI=estimated_orthobounds_CI)


#### Bounds on multidimensional effect: see 2.6.simultaneous_bounds.R

#### Bounds with logistic lasso