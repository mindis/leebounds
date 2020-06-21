rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")

library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)
library(hdm)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data_covariates<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data_covariates)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009covariates.feather"))
### Lee's covariates
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
variables_for_outcome<-baseline_varnames
Lee_data_covariates<-Lee_data[,baseline_varnames]
Lee_data_all_covariates<-cbind(Lee_data_all_covariates,Lee_data_covariates[,baseline_varnames])

source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))
weeks=1:208
selected_weeks<-c(45,90,104,135,180,208)
N_rep=1000
ci_alpha=0.05
quantile_grid_size=0.01

selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))

orthoestimates_postlasso<-matrix(0,2,length(selected_weeks))
CR_ortho_postlasso<-matrix(0,2,length(selected_weeks))

selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))
form_nonmonotone<-as.formula("selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412")

## save the results
s.hat.lasso_std<-list()
y.hat.lasso_std<-list()

form_final_selection<-list()
selected_covs_selection<-list()
### preliminary lasso run for employment equation
for (i in 1:6) {
  print(i)
  week<-selected_weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_all_covariates)
  leedata_cov[is.na(leedata_cov)]<-0
  min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
  max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  
  if (i>=4) {
    glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function = rlassologit,
                                selection_function_name="rlassologit",
                                names_to_include = baseline_varnames,standardize=FALSE)
    
  } else {
    glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function = rlassologit,selection_function_name="rlassologit",standardize=FALSE)
    
  }
  
  selected_covs_selection[[i]]<-setdiff(unique( unlist( strsplit (c(
    names(glm.fit$coefficients[glm.fit$coefficients!=0 ])
  ),split=":"))), c("(Intercept)","treat" ) )
  
  # lm.fit<-rlasso(outcome~.,le)
}
united_covs<-unique(unlist(selected_covs_selection))
selected_covariates_outcome<-list()
selected_covariates_outcome[[1]]<-c("FEMALE","BLACK","EARN_YR")
selected_covariates_outcome[[2]]<-c("FEMALE","BLACK","EARN_YR","MOSINJOB")
selected_covariates_outcome[[3]]<-c("FEMALE","BLACK","MOSINJOB")
selected_covariates_outcome[[4]]<-c("BLACK","AGE", "EARN_YR","MOSINJOB", "HH_INC5")
selected_covariates_outcome[[5]]<-c("FEMALE" , "BLACK",   "EARN_YR",   "MOSINJOB", "HH_INC5")
#selected_covariates_outcome[[5]]<-c("FEMALE" , "BLACK",   "EARN_YR","MOSINJOB" )
selected_covariates_outcome[[6]]<-c("FEMALE","BLACK")
for (i in 1:length(selected_weeks)) {
  
  print(i)
  week<-selected_weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_all_covariates)
  leedata_cov[is.na(leedata_cov)]<-0
  min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
  max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  

  glm.fit<-estimate_selection(form=form_nonmonotone,leedata=leedata_cov,selection_function = rlassologit,selection_function_name = "rlassologit")
    
  

  
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  s.hat.lasso_std[[i]]<-s.hat
  p.0.star<-s.hat.lasso_std[[i]]$s.0.hat/s.hat.lasso_std[[i]]$s.1.hat
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat.lasso_std[[i]],
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=selected_covariates_outcome[[i]],min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,p0_cuttoff=0.92) 
  
  orthoestimates_postlasso[,i]<-GetBounds( leebounds_ortho_result)

  
  y.hat<-leebounds_ortho_result$y.hat
  y.hat.lasso_std[[i]]<-y.hat
  
  set.seed(1)
  
  weights_matrix<-matrix(rexp(N_rep*dim(leedata_cov)[1]),ncol=N_rep)
  weights_list<-lapply(seq_len(ncol(weights_matrix)), function(i) weights_matrix[,i])
  estimated_bb<-lapply(weights_list,function(x) GetBounds(second_stage_wrapper(weights=x,leedata=leedata_cov,y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat,inds_helps=leebounds_ortho_result$inds_helps) ))
  estimated_bb<-matrix(unlist(estimated_bb),byrow=FALSE,nrow=2)
  
  CR_ortho_postlasso[,i]<-compute_confidence_region(ATE_boot=t(estimated_bb),
                                                    ATE_est=     orthoestimates_postlasso [,i],
                                                    ci_alpha=ci_alpha)
  
}

write.csv(orthoestimates_postlasso,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col5_estimates_single_spec.csv"))
write.csv(CR_ortho_postlasso,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col5_CR_single_spec.csv"))
