rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")

#library(sas7bdat)
#library(SDMTools)
#library(expm)

library(feather)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009covariates.feather"))


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
Lee_data_all_covariates<-cbind(Lee_data[,baseline_varnames],Lee_data_all_covariates)

weeks<-1:208
employment_rate_X1_treat<-rep(0,length(weeks))
employment_rate_X1_control<-rep(0,length(weeks))
employment_rate_X0_treat<-rep(0,length(weeks))
employment_rate_X0_control<-rep(0,length(weeks))
estimated.p.0.hat.nonmonotone<-matrix(0,dim(Lee_data_covariates),length(weeks))

selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula("selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412")

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
  
  leedata_cov_nm<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,Lee_data_all_covariates)
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

write.csv(Figure1_dataset,paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure1_selected_covs.csv"))

selected_weeks<-weeks
Figure2_dataset<-data.frame(weeks=weeks[selected_weeks],
                            fraction=fraction_treat_helps[selected_weeks])
write.csv(Figure2_dataset,paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure2_selected_covs.csv"))

#save.image(paste0(my_path,"/JobCorps/Figure_Workspaces/Figure12.RData"))

