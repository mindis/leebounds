rm(list=ls())
library(feather)
library(tidyverse)

### BSS (2019) test hyperparameter
beta.N<-0.01
## number of bootstrap repetition
B<-1000

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
Lee_data_all_covariates[is.na(Lee_data_all_covariates)]<-0
Lee_data_all_covariates$EARN_YR_quant<-Lee_data_all_covariates$EARN_YR>720 & Lee_data_all_covariates$EARN_YR<3315
Lee_data_all_covariates$EARN_YR_perc<-Lee_data_all_covariates$EARN_YR>3315 & Lee_data_all_covariates$EARN_YR<7298
source(paste0(my_path,"/JobCorps/STEP2_Test_Monotonicity/utils_for_test.R"))

week=90
B<-1000
cov_names<-c("TYPEWORR5","EARN_YR_perc")

group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=20]<-1
mygroup<-ungroup(mygroup)
mygroup$group[mygroup$group %in% c(0,2,3)]<-0

pvalue<-test_wrapper(week,mygroup)

