rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")

#library(sas7bdat)
#library(SDMTools)
#library(expm)


my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)


#Lee_data<-read_feather(paste0(my_path,"/R_code/dataLee2009.feather"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))

TREAT<-Lee_data$TREATMNT.y
# covariates selected by David Lee (2008)
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

weeks<-1:208
employment_rate_X1_treat<-rep(0,length(weeks))
employment_rate_X1_control<-rep(0,length(weeks))
employment_rate_X0_treat<-rep(0,length(weeks))
employment_rate_X0_control<-rep(0,length(weeks))
estimated.p.0.hat.nonmonotone<-matrix(0,dim(Lee_data_covariates),length(weeks))

selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula(paste0("selection~treat*(", paste0(selected_names,collapse="+"),")"))

weeks<-1:208
s.hat.nonmonotone.glm<-array(0,c(dim(Lee_data_covariates)[1],2,length(weeks)))
fraction_treat_helps<-rep(0,length(weeks))
for (i in 1:length(weeks)) {
  week<-weeks[i]
  print(paste0("Week ",i))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_cov_nm<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,Lee_data_covariates)
  leedata_cov_nm[is.na(leedata_cov_nm)]<-0
  
  glm.fit<-glm(form=form_nonmonotone,data=leedata_cov_nm,family="binomial")
  res<-predict_selection(glm.fit,leedata_cov_nm)
  s.hat = data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  
  
  estimated.p.0.hat.nonmonotone[,i]<-s.hat$s.0.hat/s.hat$s.1.hat
  
  employment_rate_X0_treat[i]<-mean(leedata_cov_nm$selection[leedata_cov_nm$treat==1 & estimated.p.0.hat.nonmonotone[,i]<1 ])
  employment_rate_X0_control[i]<-mean(leedata_cov_nm$selection[leedata_cov_nm$treat==0 & estimated.p.0.hat.nonmonotone[,i]<1 ])
  
  employment_rate_X1_treat[i]<-mean(leedata_cov_nm$selection[leedata_cov_nm$treat==1 & estimated.p.0.hat.nonmonotone[,i]>1 ])
  employment_rate_X1_control[i]<-mean(leedata_cov_nm$selection[leedata_cov_nm$treat==0 & estimated.p.0.hat.nonmonotone[,i]>1 ])
  s.hat.nonmonotone.glm[,,i]<-as.matrix(s.hat,ncol=2)
  fraction_treat_helps[i]<-mean(estimated.p.0.hat.nonmonotone[,i]<=1)
}



selected_inds<-weeks
Figure1_dataset=data.frame(weeks=c(weeks[selected_inds],weeks[selected_inds]),
                           delta=c((employment_rate_X0_treat-employment_rate_X0_control)[selected_inds],
                                   (employment_rate_X1_treat-employment_rate_X1_control)[selected_inds]),
                           group=c(rep(1,length(selected_inds)),rep(0,length(selected_inds)))
)

write.csv(Figure1_dataset,paste0(my_path,"/JobCorps/STEP3_MISC/csv/Figure1.csv"))

selected_weeks<-weeks
Figure2_dataset<-data.frame(weeks=weeks[selected_weeks],
                            fraction=fraction_treat_helps[selected_weeks])
write.csv(Figure2_dataset,paste0(my_path,"/JobCorps/STEP3_MISC/csv/Figure2.csv"))

#save.image(paste0(my_path,"/JobCorps/Figure_Workspaces/Figure12.RData"))

