## Step 3.

## I assume that monotonicity holds S_{0i} >= S_{1i} (treatment can only hurt for employment)
## for weeks 20-80 where the trimming probability is bounded away from zero

## I report Lee bounds with covariates, 95% pointwise bounds for the TE, and 95% uniform bounds for the treatment effect function 
## of weeks.

rm(list=ls())

args<-commandArgs(TRUE)

if (length(args)<1) {
  weeks = 1:208
  selection_function_name="glm"
  selection_function=glm
  
} else {
  min_week = as.numeric(args[1])
  max_week = as.numeric(args[2])
  weeks = min_week:max_week
  selection_function_name=args[3]
  
  
}

#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
#install.packages("feather",repo="http://cran.rstudio.com/")
#library(feather)
library(quantreg)
library(expm)
library(purrr)


#install.packages("rlist")
#library(rlist)

source(paste0(my_path,"/R/leebounds.R"))
#source(paste0(my_path,"/R/libraries.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
if (selection_function_name=="rlassologit") {
  Lee_data_covariates<-read_feather(paste0(my_path,"/R/dataLee2009covariates.feather"))
  library(hdm)
} else {
  ## same covariates as in Lee (2009)
  baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                       "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                       "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                       "PERS_INC2","PERS_INC3","PERS_INC4",
                       
                       "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                       "WKEARNR")
  
  Lee_data_covariates<-Lee_data[,baseline_varnames]
}

sample_size<-dim(Lee_data_covariates)[1]
form_monotone<-"selection~."
selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula(paste0("selection~treat*(", paste0(selected_names,collapse="+"),")"))
selection_function=get(selection_function_name)


### Initialize arrays #####
s.0.hat.monotone<-matrix(0,sample_size,length(weeks))
s.0.hat.nonmonotone<-matrix(0,sample_size,length(weeks))
s.1.hat.monotone<-matrix(0,sample_size,length(weeks))
s.1.hat.nonmonotone<-matrix(0,sample_size,length(weeks))
p.0.hat.monotone<-matrix(0,sample_size,length(weeks))
p.0.hat.nonmonotone<-matrix(0,sample_size,length(weeks))

colnames(s.0.hat.monotone)<-weeks
colnames(s.1.hat.monotone)<-weeks
colnames(p.0.hat.monotone)<-weeks

colnames(s.0.hat.nonmonotone)<-weeks
colnames(s.1.hat.nonmonotone)<-weeks
colnames(p.0.hat.nonmonotone)<-weeks

lasso_res<-list()
lasso_res_nm<-list()

for (j in 1:length(weeks)) {
  week<-weeks[j]
  print(paste0("Estimating selection equation for week ",week))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leedata_cov<-cbind(leedata_week,Lee_data_covariates)
  leedata_cov[is.na(leedata_cov)]<-0
  
  lasso.fit<-estimate_selection(leedata_cov,selection_function=selection_function,form=form_monotone)
  res<-predict_selection(lasso.fit,leedata_cov)
  lasso_res[[j]]<-lasso.fit$coefficients
  
  
  
  # lasso_res[[week]]<-sapply(lasso.fit$coefficients,round,4)
  
  s.0.hat.monotone[,j]<-res$s.0.hat
  s.1.hat.monotone[,j]<-res$s.1.hat
  p.0.hat.monotone[,j]<-s.0.hat.monotone[,j]/s.1.hat.monotone[,j]
  
  lasso.fit.nm<-estimate_selection(leedata_cov,selection_function=selection_function,form= form_nonmonotone)
  lasso_res_nm[[j]]<-lasso.fit.nm$coefficients
  
  res_nm<-predict_selection( lasso.fit.nm,leedata_cov)
  s.0.hat.nonmonotone[,j]<- res_nm$s.0.hat
  s.1.hat.nonmonotone[,j]<- res_nm$s.1.hat
  p.0.hat.nonmonotone[,j]<-s.0.hat.nonmonotone[,j]/s.1.hat.nonmonotone[,j]
  
  # lasso_res_nm[[week]]<-sapply(lasso.fit.nm$coefficients,round,4)
}


write.csv(as.data.frame(s.0.hat.monotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/s.0.hat.monotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))
write.csv(as.data.frame(s.1.hat.monotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/s.1.hat.monotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))

write.csv(as.data.frame(s.0.hat.nonmonotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/s.0.hat.nonmonotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))
write.csv(as.data.frame(s.1.hat.nonmonotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/s.1.hat.nonmonotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))

write.csv(as.data.frame(p.0.hat.nonmonotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/p.0.hat.nonmonotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))
write.csv(as.data.frame(p.0.hat.monotone),paste0("First_Stage_Predicted_Values/Predicted_Selection/p.0.hat.monotone.",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".csv"))





rm(list=c("Lee_data","leedata_cov","leedata_week","Lee_data_covariates"))
setwd(paste0(my_path,"/R/2.estimate_bounds/"))
save.image(paste0("First_Stage_Predicted_Values/estimated_selection_",selection_function_name,"_weeks_",min(weeks),"_",max(weeks),".RData"))

#savehistory(file="rlassologit.RData")