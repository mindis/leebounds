rm(list=ls())

args<-commandArgs(TRUE)

if (length(args)<1) {
  # default choice: all weeks, 0.01 grid size for quantile regression
  weeks = 1:208
  quantile_grid_size=0.01
} else {
  min_week = as.numeric(args[1])
  max_week = as.numeric(args[2])
  weeks = min_week:max_week
  quantile_grid_size=as.numeric(args[3])
}
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"~"
library(feather)
library(quantreg)
library(expm)
library(purrr)
library(hdm)

my_path<-"/net/holyparkesec/data/tata/Vira/"
source(paste0(my_path,"/R_code/leebounds.R"))
#source(paste0(my_path,"/R_code/libraries.R"))
source(paste0(my_path,"/R_code/ortholeebounds.R"))
Lee_data<-read_feather(paste0(my_path,"/R_code/dataLee2009.feather"))
#Lee_data_covariates<-read_feather(paste0(my_path,"/R_code/dataLee2009covariates.feather"))
setwd(paste0(my_path,"/R_code/2.estimate_bounds/"))

sample_size<-dim(Lee_data)[1]
taus<-seq(0,1,by=quantile_grid_size)
estimated_quantiles_11=array(0,c(sample_size,length(taus),length(weeks)))
estimated_quantiles_10=array(0,c(sample_size,length(taus),length(weeks)))

  
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")

leedata_week<-list()
for (j in 1:length(weeks)) {
  week<-weeks[j]
  print(week)
  print(paste0("Estimating quantile table for week ",week))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week[[j]]<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind( leedata_week[[j]],Lee_data[,baseline_varnames])
  leedata_cov[is.na(leedata_cov)]<-0
  

  
    estimated_quantiles_11[,,j]<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==1 & leedata_cov$selection==1,],
                                                                  test_data=leedata_cov,
                                                                  variables_for_outcome=baseline_varnames,
                                                                  taus=taus)
    estimated_quantiles_10[,,j]<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                                                  test_data=leedata_cov,
                                                                  variables_for_outcome=baseline_varnames,
                                                                  taus=taus)
    
  
 

}
rm(list=c("Lee_data","leedata_cov"))
save.image(paste0("First_Stage_Predicted_Values/estimated_quantiles_rq_",as.character(quantile_grid_size),"_weeks_",min(weeks),"_",max(weeks),".RData"))
