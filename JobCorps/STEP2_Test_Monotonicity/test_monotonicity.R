rm(list=ls())
library(feather)
library(tidyverse)
### for parallel computing, replacable by standard for loop
library(doParallel)
library(foreach)
cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))

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



myres_M0=foreach(week=60:79, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  
  
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
  mygroup$group[mygroup$group %in% c(0,1,2)]<-0
  
  pvalue<-test_wrapper(week,mygroup)
  
}


myres_M1=foreach(week=80:108, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  
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
  
}


myres_M2=foreach(week=109:139, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
 
  if (week<=111) {
    cov_names<-c("REASED_R4","WELF5AC")
    
    group_weight<-1:length(cov_names)
    mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
    mygroup<-as.numeric(mygroup)
    
    mygroup<-group_by(data.frame(group=mygroup),group) %>%
      count %>%
      inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
    ## unite groups of size less than 30 into 1 group
    mygroup$group[mygroup$n<=20]<-1
    mygroup<-ungroup(mygroup)
  } else {
    
 
  cov_names<-c("EARN_YR_quant","R_HOME1","WELF5AC","TYPEWORR5")
 
  group_weight<-1:length(cov_names)
  mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
  mygroup<-as.numeric(mygroup)
  
  mygroup<-group_by(data.frame(group=mygroup),group) %>%
    count %>%
    inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
  ## unite groups of size less than 30 into 1 group
  mygroup$group[mygroup$n<=20]<-1
  mygroup<-ungroup(mygroup)
  mygroup$group[mygroup$group %in% c(0,1,2,7)]<-0
  mygroup$group[mygroup$group %in% c(3,6)]<-3
  }
  pvalue<-test_wrapper(week,mygroup)
}

myres_M3=foreach(week=140:153, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  cov_names<-c("REASED_R4","IMP_PRO1")
  
  
  group_weight<-1:length(cov_names)
  mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
  mygroup<-as.numeric(mygroup)
  
  mygroup<-group_by(data.frame(group=mygroup),group) %>%
    count %>%
    inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
  ## unite groups of size less than 30 into 1 group
  mygroup$group[mygroup$n<=20]<-1
  mygroup<-ungroup(mygroup)
 
  
  pvalue<-test_wrapper(week,mygroup)
  
}



myres_M4=foreach(week=154:164, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  cov_names<-c("REASED_R4","IMP_PRO1")
  
  
  group_weight<-1:length(cov_names)
  mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
  mygroup<-as.numeric(mygroup)
  
  
  mygroup<-group_by(data.frame(group=mygroup),group) %>%
    count %>%
    inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
  ## unite groups of size less than 30 into 1 group
  mygroup$group[mygroup$n<=20]<-1
  mygroup<-ungroup(mygroup)
  
  mygroup$group[mygroup$group %in% c(0,2)]<-0
  if (week <=156) {
    mygroup$group[mygroup$group %in% c(1,3)]<-1
  }

  pvalue<-test_wrapper(week,mygroup)
  
}



myres_M5=foreach(week=165:180, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  cov_names<-c("FRQ_POT3","BLACK")
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
  
}

myres_M6=foreach(week=181:208, .combine = 'cbind',.packages = c("tidyverse"))  %dopar% {
  ## constructed data for week
  # cov_names<-c("REASED_R4","IMP_PRO1")
  cov_names<-c("DRG_SUMP2","IMPRCNTR4","TYPEJOBB9")
  group_weight<-1:length(cov_names)
  mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
  mygroup<-as.numeric(mygroup)
  
  mygroup<-group_by(data.frame(group=mygroup),group) %>%
    count %>%
    inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID))
  ## unite groups of size less than 30 into 1 group
  mygroup$group[mygroup$n<=20]<-1
  mygroup<-ungroup(mygroup)
  mygroup$group[mygroup$group %in% c(0,1,2,3)]<-0
  
  pvalue<-test_wrapper(week,mygroup)
  
}

test_result<-data.frame(weeks=c(1:208),pvalue=rep(1,208))
test_result$pvalue[60:79]<-as.numeric(myres_M0)
test_result$pvalue[80:108]<-as.numeric(myres_M1)
test_result$pvalue[109:139]<-as.numeric(myres_M2)
test_result$pvalue[140:153]<-as.numeric(myres_M3)
test_result$pvalue[154:164]<-as.numeric(myres_M4)
test_result$pvalue[165:180]<-as.numeric(myres_M5)
test_result$pvalue[181:208]<-as.numeric(myres_M6)

write.csv(test_result,paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/test_result.csv"))