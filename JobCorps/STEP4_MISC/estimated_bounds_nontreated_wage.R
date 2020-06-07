rm(list=ls())
#install.packages("feather")
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
source(paste0(my_path,"/R/orthogonal_correction.R"))
### FIGURE 1 and 2###
weeks<-1:208
quantile_grid_size=0.01
### FIGURE 5
s.hat.nonmonotone.glm<-array(0,c(dim(Lee_data_all_covariates)[1],2,length(weeks)))
estimated.p.0.hat.nonmonotone<-matrix(0,dim(Lee_data_all_covariates)[1],length(weeks))
bounds_non_treated_wage<-matrix(0,2,length(weeks))
frac_positive<-rep(0,length(weeks))
estimates_plb<-matrix(0,2,length(weeks))



TREAT<-Lee_data$TREATMNT.y
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR"  )

Lee_data_all_covariates<-cbind(Lee_data[,baseline_varnames],Lee_data_all_covariates)


selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula("selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412")
#form_nonmonotone<-as.formula(paste0("selection~(treat)*(", paste0(baseline_varnames,collapse="+"),")"))

s.hat.nonmonotone.glm<-array(0,c(dim(Lee_data_all_covariates)[1],2,length(weeks)))
fraction_treat_helps<-rep(0,length(weeks))
variables_for_outcome<-baseline_varnames




for (i in 1:length(weeks)) {
  print(i)
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_all_covariates)
  leedata_cov[is.na( leedata_cov)]<-0
  ## analysis without monotonicity
  glm.fit<-glm(form=form_nonmonotone,data= leedata_cov[,setdiff(colnames(leedata_cov),"outcome")],family="binomial")
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  s.hat.nonmonotone.glm[,,i]<-as.matrix(s.hat,ncol=2)
  estimated.p.0.hat.nonmonotone[,i]<-s.hat$s.0.hat/s.hat$s.1.hat
  
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat,
                                          variables_for_outcome=baseline_varnames,
                                          quantile_grid_size=quantile_grid_size,
                                          bounds_fun=ortho_bounds_nontreated_wage_ss,weights=Lee_data$DSGN_WGT.y,c_quant=1,ortho=TRUE)
  
  bounds_non_treated_wage[,i]<-GetBounds(leebounds_ortho_result)
  #nonortho_bounds[,i]<- leebounds_ortho_result$
}

#selected_weeks<-4:208
#selected_weeks<-5:208

selected_weeks<-5:208
Figure5_dataset<-data.frame(weeks=c(weeks[selected_weeks],
                                    weeks[selected_weeks]),
                            bound=c(bounds_non_treated_wage[1,selected_weeks],bounds_non_treated_wage[2,selected_weeks]),
                            group=c(rep(0,length(selected_weeks)),rep(1,length(selected_weeks))))
write.csv(Figure5_dataset,paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure5.csv"))