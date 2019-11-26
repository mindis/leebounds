# Table 4
rm(list=ls())
library(sas7bdat)
library(dplyr)
library(tidyverse)
library(SDMTools)
library(feather)
library(expm)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/Data/dataLee2009.csv"))


Lee_data<-as.data.frame(Lee_data)

#Lee_data<-read_feather(paste0(my_path,"/R_code/dataLee2009.feather"))
source(paste0(my_path,"/leebounds.R"))
source(paste0(my_path,"/utils.R"))

TREAT<-Lee_data$TREATMNT.y
EMPLOYED208<-Lee_data$HWH208>0
logWAGE208<-log(Lee_data$EARNH208/Lee_data$HWH208)
logWAGE208[is.na(logWAGE208)]<-0
logWAGE208[logWAGE208==-Inf]<-0
leedata<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=EMPLOYED208,outcome=logWAGE208))
int_weights<-round(Lee_data$DSGN_WGT.y)
weighted_leedata<-data.frame()
for (wt in unique(int_weights)) {
  subset<-leedata[int_weights==wt,]
  subset_wt<-subset
  for (i in 1:(wt-1)){
    subset_wt<-rbind(subset_wt,subset)
  }
  weighted_leedata<-rbind(weighted_leedata,subset_wt)
}
leebounds_result_weighted<-leebounds(weighted_leedata)
JobCorps_bounds_weighted<-GetBounds(leebounds_result_weighted)

leebounds_result<-leebounds_unknown_sign(leedata)
JobCorps_bounds_unweighted<-GetBounds(leebounds_result)
### bootstrap standard errors

JobCorps_bb<-main_bb(leedata,N_rep=1000,function_name=leebounds_unknown_sign)
JobCorps_se<-apply(JobCorps_bb,2,sd)
JobCorps_ci<-c(JobCorps_bounds_unweighted[1]-1.96*JobCorps_se[1],
               JobCorps_bounds_unweighted[2]+1.96*JobCorps_se[2])
JobCorps_ci_ImbensManski<-c(JobCorps_bounds_unweighted[1]-1.645*JobCorps_se[1],
                            JobCorps_bounds_unweighted[2]+1.645*JobCorps_se[2])
JobCorps_optimal_ci<-compute_confidence_region(JobCorps_bb,JobCorps_bounds_unweighted,ci_alpha=0.05)

table3<-cbind(JobCorps_bounds_unweighted,JobCorps_se,JobCorps_ci,JobCorps_ci_ImbensManski,JobCorps_optimal_ci)
table3<-apply(table3,2,round,3)
write.csv(table3,paste0(my_path,"/1.replicate_lee2009/tables/Table3_replicated.csv"))
## Table 6
weeks<-c(45,90,135,180)
JobCorps_bounds_unweighted_table6<-matrix(0,nrow=2,ncol=length(weeks))
for (i in 1:length(weeks)) {
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  leebounds_week<-leebounds_unknown_sign(leedata_week)
  JobCorps_bounds_unweighted_table6[,i]<-GetBounds( leebounds_week)

}
JobCorps_bounds_unweighted_table6<-apply(JobCorps_bounds_unweighted_table6,2,round,3)
write.csv(JobCorps_bounds_unweighted_table6,paste0(my_path,"/1.replicate_lee2009/tables/Table6_replicated.csv"))
