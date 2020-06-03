rm(list=ls())

## load packages
print ("Loading packages ...")
library(reldist)
library(expm)
library(feather)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
sink(paste0(my_path,"/JobCorps/STEP2_Test_Monotonicity/Table_Replicate_Lee.log"))
print ("Loading data ...")
Lee_data_covariates<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data_covariates)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps_data/dataLee2009covariates.feather"))

baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_all_covariates<-cbind(Lee_data_all_covariates,Lee_data[,baseline_varnames])
Lee_data_all_covariates[is.na(Lee_data_all_covariates)]<-0

Lee_data_all_covariates$EARN_YR_perc<-Lee_data_all_covariates$EARN_YR>3315 & Lee_data_all_covariates$EARN_YR<7298

source(paste0(my_path,"/JobCorps/STEP2_Test_Monotonicity/utils_for_test.R"))
source(paste0(my_path,"/JobCorps/STEP2_Test_Monotonicity/utils_for_bounds.R"))


selected_weeks<-c(45,90,104,135,180,208)
N_rep=1000
ci_alpha=0.05
quantile_grid_size=0.01
# Column 1
estimates<-matrix(0,2,length(selected_weeks))
CR<-matrix(0,2,length(selected_weeks))
IM_CR<-matrix(0,2,length(selected_weeks))

# Column 2
estimates_sharp<-matrix(0,2,length(selected_weeks))
CR_sharp<-matrix(0,2,length(selected_weeks))
IM_CR_sharp<-matrix(0,2,length(selected_weeks))

# Column 3
estimates_sharp2<-matrix(0,2,length(selected_weeks))
CR_sharp2<-matrix(0,2,length(selected_weeks))
IM_CR_sharp2<-matrix(0,2,length(selected_weeks))


## week 208 predicted wage potential

logwage_week<-as.numeric(as.matrix(log(Lee_data$EARNH208/Lee_data$HWH208)))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0
leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_covariates)
leedata_cov[is.na( leedata_cov)]<-0

lm.fit<-lm(form=as.formula(paste0("outcome~",paste0(baseline_varnames,collapse="+")) ), 
           data=leedata_cov[leedata_cov$selection==1,])
Lee_data$predicted_log_wage<-predict(lm.fit,leedata_cov)
Lee_data$predicted_wage_group<-sapply(exp(Lee_data$predicted_log_wage),Lee_grouping)


## week 208 predicted ML wage potential
logwage_week<-as.numeric(as.matrix(log(Lee_data$EARNH208/Lee_data$HWH208)))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0
leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_all_covariates)
leedata_cov[is.na( leedata_cov)]<-0

selected_names<-setdiff(colnames(Lee_data_all_covariates),"MPRID")
glm.fit<-rlassologit(form=as.formula(paste0("selection~",paste0(selected_names,collapse="+"))) , 
           data=leedata_cov)
Lee_data$predicted_emp<-predict(glm.fit,leedata_cov,type="response")
Lee_data$predicted_wage_group_ml<-sapply(Lee_data$predicted_emp,Lee_grouping,cutoffs=c(0.50,0.60,0.70,0.80))
Lee_data$predicted_wage_group_ml[is.na(Lee_data$predicted_wage_group_ml)]<-4

for (i in 1:length(selected_weeks)) {
  # prepare data
  week<-selected_weeks[i]
  print (paste0("Results for week ", week))
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  
  leedata_week<-Lee_data[,c("MPRID","TREATMNT.y","DSGN_WGT.y",earn_name,hwh_name)] %>%
    inner_join (Lee_data[,c("MPRID","predicted_wage_group")],by= c("MPRID"="MPRID")) %>%
    inner_join (Lee_data[,c("MPRID","predicted_wage_group_ml")],by= c("MPRID"="MPRID"))
  
  colnames(leedata_week)<-c("MPRID","treat","weights",earn_name,hwh_name,"predicted_wage_group","predicted_wage_group_ml")
  leedata_week$selection<-leedata_week[,earn_name]>0
  leedata_week$outcome<-as.numeric(as.matrix(log(leedata_week[,earn_name]/leedata_week[,hwh_name])))
  ### COLUMN 1
  
  print ("Computing  Lee (2009) estimates (column 1) ...")
  ## compute bounds 
  estimates[,i]<-GetBounds(leebounds_unknown_sign(leedata_week))
  ##  bootstrap draws of bounds
  print ("Computing  Lee (2009) confidence region by regular bootstrap (column 1) ...")
  
 # bounds_bb<-main_bb(function_name=leebounds_unknown_sign,mydata=leedata_week,N_rep=N_rep)
  # ## confidence region for identified set
  #CR[,i]<-compute_confidence_region(bounds_bb,estimates[,i], ci_alpha=ci_alpha )
  #IM_CR[,i]<-imbens_manski(bounds_bb,estimates[,i], ci_alpha=ci_alpha)
  
  ### Column 2 
  print ("Computing  Lee (2009) estimates (column 2) ...")
  res<-Lee_sharp_bounds(leedata_week,group_name="predicted_wage_group",treat_helps = (week>=90))
  estimates_sharp[,i]<-GetBounds(res)
  
  res2<-Lee_sharp_bounds(leedata_week,group_name="predicted_wage_group_ml",treat_helps = (week>=90))
  estimates_sharp2[,i]<-GetBounds(res2)
  
  
 # bounds_bb<-main_bb(function_name=Lee_sharp_bounds,mydata=leedata_week,N_rep=N_rep,group_name="predicted_wage_group",treat_helps = (week>=90))
  ## confidence region for identified set
#  CR_sharp[,i]<-compute_confidence_region(bounds_bb,estimates_sharp[,i], ci_alpha=ci_alpha )
 # IM_CR_sharp[,i]<-imbens_manski(bounds_bb,estimates_sharp[,i], ci_alpha=ci_alpha)
  
}

estimates_table<-rbind(estimates,estimates_sharp)
estimates_table<-t(estimates_table)
colnames(estimates_table)<-c("Lee_2009_lb","Lee_2009_ub","Lee_2009_lb_cov","Lee_2009_ub_cov")
CR_table<-rbind(CR,CR_sharp)
CR_table<-t(CR_table)

IM_CR_table<-rbind(IM_CR,IM_CR_sharp)
IM_CR_table<-t(IM_CR_table)


estimates_table<-apply(estimates_table,2,round,3)
CR_table<-apply(CR_table,2,round,3)
IM_CR_table<-apply(IM_CR_table,2,round,3)

print("Saving estimates in STEP2_Estimate_Bounds/csv/ ...")
sink(file=NULL)
closeAllConnections()
write.csv(estimates_table,paste0(my_path,"JobCorps/STEP2_Test_Monotonicity/csv/Lee_replicated_estimates.csv"))
write.csv(CR_table,paste0(my_path,"JobCorps/STEP2_Test_Monotonicity/csv/Lee_replicated_CR.csv"))
write.csv(IM_CR_table,paste0(my_path,"JobCorps/STEP2_Test_Monotonicity/csv/Lee_replicated_CR_IM.csv"))
