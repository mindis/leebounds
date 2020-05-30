rm(list=ls())
library(feather)
library(tidyverse)
### for parallel computing, replacable by standard for loop
library(doParallel)
library(foreach)
cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))

my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009covariates.feather"))

baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_all_covariates<-cbind(Lee_data_all_covariates,Lee_data[,baseline_varnames])
Lee_data_all_covariates$EARN_YR_quant<-Lee_data_all_covariates$EARN_YR>720 & Lee_data_all_covariates$EARN_YR<3315

source(paste0(my_path,"/JobCorps/STEP3_MISC/utils_for_test.R"))

### create covariate groups


cov_names<-c("PERS_INC3","EARN_YR_quant","MOS_AFDC8")
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)
mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=20]<-2
mygroup<-ungroup(mygroup)

### BSS (2019) test hyperparameter
beta.N<-0.01
## number of bootstrap repetition
B<-500



myres_M0=foreach(week=60:89, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_cov<-data.frame(treat=Lee_data$TREATMNT.y,
                          selection = logwage_week>0,
                          MPRID=Lee_data$MPRID) %>%
    inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
  
  
  
  leedata_cov$group_copy<-leedata_cov$group
  ## created grouped data
  grouped_data<-group_by(leedata_cov,group)
  grouped_data<-do.call(rbind,group_map(grouped_data,compute_mean_sd))
  grouped_data$group<-grouped_data$group_copy
  
  
  ### compute test statistic
  grouped_data<-group_by(grouped_data,group)
  grouped_data_tstat<-grouped_data
  grouped_data_tstat$MEAN<-0
  tstat<-unlist(group_map(grouped_data_tstat,compute_tstat_by_group,sign=1,mean_name="MEAN"))
  max.t.stat<-max(max(tstat),0)
  

  ## compute critical value
  
  res<-compute_crit_val2(grouped_data,tstat= max.t.stat,B=B,beta.N=beta.N)
  pvalue<- res$pvalue

}



cov_names<-c("EARN_YR_quant","R_HOME1","WELF5AC","TYPEWORR5")

set.seed(1)
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=20]<-2
mygroup<-ungroup(mygroup)
### testing M1 (Treatment helps) for weeks after 90

myres_M1=foreach(week=90:139, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_cov<-data.frame(treat=1-Lee_data$TREATMNT.y,
                          selection = logwage_week>0,
                          MPRID=Lee_data$MPRID) %>%
    inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
  
  
  
  
  leedata_cov$group_copy<-leedata_cov$group
  ## created grouped data
  grouped_data<-group_by(leedata_cov,group)
  grouped_data<-do.call(rbind,group_map(grouped_data,compute_mean_sd))
  grouped_data$group<-grouped_data$group_copy
  
 
  ### compute test statistic
  grouped_data<-group_by(grouped_data,group)
  grouped_data_tstat<-grouped_data
  grouped_data_tstat$MEAN<-0
  tstat<-unlist(group_map(grouped_data_tstat,compute_tstat_by_group,sign=1,mean_name="MEAN"))
  max.t.stat<-max(max(tstat),0)
  
  res<-compute_crit_val2(grouped_data,tstat= max.t.stat,B=B,beta.N=beta.N)
  pvalue<- res$pvalue

}
### save results
test_result<-data.frame(weeks=c(1:208),pvalue=rep(1,208))
test_result$pvalue[60:89]<-as.numeric(myres_M0)
test_result$pvalue[90:139]<-as.numeric(myres_M1)

write.csv(test_result,paste0(my_path,"/JobCorps/STEP3_MISC/test_result.csv"))