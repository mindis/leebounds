rm(list=ls())

#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"~"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/leebounds.R"))
source(paste0(my_path,"/libraries.R"))
Lee_data<-read.csv(paste0(my_path,"/Data/dataLee2009.csv"))
setwd(paste0(my_path,"/2.estimate_bounds/"))

############ WEIGHTED DATA ###############################
int_weights<-round(Lee_data$DSGN_WGT.y)
weighted_leedata<-data.frame()
for (wt in unique(int_weights)) {
  subset<-Lee_data[int_weights==wt,]
  subset_wt<-subset
  for (i in 1:(wt-1)){
    subset_wt<-rbind(subset_wt,subset)
  }
  weighted_leedata<-rbind(weighted_leedata,subset_wt)
}

baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")

# Step 1. Plot trimming proportion against week since random assignment
# all weeks
weeks<-c(1:208)
# weeks from Lee Bounds
selected_weeks<-c(45,90,135,180)



prop_trim<-rep(0,length(weeks))
prop_trim_weighted<-rep(0,length(weeks))
estimated_weighted_leebounds<-matrix(0,2,length(weeks))
estimated_leebounds<-matrix(0,2,length(weeks))

for (i in 1:length(weeks)) {
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  logwage_week[logwage_week>0]<-logwage_week[logwage_week>0]+rnorm(sum(logwage_week>0),0.001)
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind(leedata_week,Lee_data[,baseline_varnames])
  ## treating NAs in baseline covariates
  leedata_cov[is.na(leedata_cov)]<-0
  print(i)
  prop_trim[i]<-mean(leedata_cov$selection==1 & leedata_cov$treat==0)/mean(leedata_cov$selection==1 & leedata_cov$treat==1)
  
}


png("Figures/Trimming_unweighted.png")
plot(weeks,prop_trim,xlab="Weeks since random assignment",ylab="Trimming threshold",lwd=3,ylim=c(0,2),type="l")+
  abline(a=0,b=0,col="black",lwd=2)+
  title("Trimming threshold")+
  points(x=selected_weeks,y=prop_trim[selected_weeks],col="black",lwd=10)
dev.off()
################## weighted data ###############################



for (i in 1:length(weeks)) {
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(weighted_leedata[,earn_name]/weighted_leedata[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  logwage_week[logwage_week>0]<-logwage_week[logwage_week>0]+rnorm(sum(logwage_week>0),0.001)
  
  leedata_week<-as.data.frame(cbind(treat=weighted_leedata$TREATMNT.y,selection=as.numeric(weighted_leedata[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind(leedata_week,weighted_leedata[,baseline_varnames])
  ## treating NAs in baseline covariates
  leedata_cov[is.na(leedata_cov)]<-0
  print(i)
  prop_trim_weighted[i]<-mean(leedata_cov$selection==1 & leedata_cov$treat==0)/mean(leedata_cov$selection==1 & leedata_cov$treat==1)
  
}

png("Figures/Trimming_weighted.png")
plot(weeks,prop_trim_weighted,xlab="Weeks since random assignment",ylab="Trimming threshold",lwd=3,type="l",ylim=c(0.5,2))+
  abline(a=0,b=0,col="black",lwd=2)+
  title("Trimming threshold")+
  points(x=selected_weeks,y=prop_trim_weighted[selected_weeks],col="black",lwd=10)
dev.off()
####### Step 2. Report unconditional Lee bounds for each week



for (i in 1:length(weeks)) {
  week<-weeks[i]
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==0))
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==1))
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==0)/mean(leedata_cov$selection==1 & leedata_cov$treat==1))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  logwage_week[logwage_week>0]<-logwage_week[logwage_week>0]+rnorm(sum(logwage_week>0),0.001)
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind(leedata_week,Lee_data[,baseline_varnames])
  ## treating NAs in baseline covariates
  leedata_cov[is.na(leedata_cov)]<-0
  leebounds<-leebounds_unknown_sign(leedata_week)
  estimated_leebounds[,i]<-GetBounds(leebounds)
}
png("Figures/Basic_Lee_bounds.png")
plot(weeks,estimated_leebounds[1,],xlab="Weeks since random assignment",ylab="Lower and upper bound",col="blue",
     ylim=c(-1,1),type="l",lwd=3)+
  lines(estimated_leebounds[2,],col="blue",lwd=3)+
  title("Lee (2009) bounds")+
  points(x=selected_weeks,y=estimated_leebounds[1,selected_weeks],col="black",lwd=10)+
  points(x=selected_weeks,y=estimated_leebounds[2,selected_weeks],col="black",lwd=10)+
  abline(a=0,b=0,col="black",lwd=2)
dev.off()




####### Step 2. Report weighted unconditional Lee bounds for each week (weighted)



for (i in 1:length(weeks)) {
  week<-weeks[i]
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==0))
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==1))
  print(mean(leedata_cov$selection==1 & leedata_cov$treat==0)/mean(leedata_cov$selection==1 & leedata_cov$treat==1))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(weighted_leedata[,earn_name]/weighted_leedata[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  logwage_week[logwage_week>0]<-logwage_week[logwage_week>0]+rnorm(sum(logwage_week>0),0.001)
  
  leedata_week<-as.data.frame(cbind(treat=weighted_leedata$TREATMNT.y,selection=as.numeric(weighted_leedata[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind(leedata_week,weighted_leedata[,baseline_varnames])
  ## treating NAs in baseline covariates
  leedata_cov[is.na(leedata_cov)]<-0
  leebounds_weighted<-leebounds_unknown_sign(leedata_week)
  estimated_weighted_leebounds[,i]<-GetBounds(leebounds_weighted)
}
png("Figures/Basic_weighted_Lee_bounds.png")
plot(weeks, estimated_weighted_leebounds[1,],xlab="Weeks since random assignment",ylab="Lower and upper bound",col="blue",type="l",
     ylim=c(-1,1),lwd=3)+
  lines( estimated_weighted_leebounds[2,],col="blue",lwd=3)+
  abline(a=0,b=0,col="black",lwd=2)+
  title("Lee (2009) bounds")+
  points(x=selected_weeks,y= estimated_weighted_leebounds[1,selected_weeks],col="black",lwd=3)+
  points(x=selected_weeks,y= estimated_weighted_leebounds[2,selected_weeks],col="black",lwd=3)
  
dev.off()

  
