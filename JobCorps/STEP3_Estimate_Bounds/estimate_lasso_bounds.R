rm(list=ls())

## load packages
print ("Loading packages ...")
library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)
library(reldist)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
#sink(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/Table1_Col123.log"))
print ("Loading data ...")
Lee_data_covariates<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps_data/dataLee2009covariates.feather"))

Lee_data<-as.data.frame(Lee_data_covariates)

source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/orthogonal_correction.R"))
source(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/utils.R"))
selected_weeks<-c(45,90,104,135,180,208)
N_rep=1000
ci_alpha=0.05
quantile_grid_size=0.01

orthoestimates_postlasso<-matrix(0,2,length(selected_weeks))
CR_ortho<-matrix(0,2,length(selected_weeks))
IM_ortho<-matrix(0,2,length(selected_weeks))

estimates_plb<-matrix(0,2,length(selected_weeks))
frac_positive<-matrix(0,2,length(selected_weeks))

#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

z="selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412"
selected_covs<-setdiff(unlist(strsplit(z,"+",fixed=TRUE)),"selection~treat")
selected_covs_treat<-grep("treat:",selected_covs,value=TRUE)
selected_names<-unique(setdiff(c( setdiff(selected_covs,selected_covs_treat),
                                  colnames(Lee_data_covariates)),c("treat","selection","outcome","(Intercept)","X.Intercept.")))
form_nonmonotone_ss<-as.formula(paste0("selection~treat+(",paste0(selected_names,collapse="+"), ")*(",paste0(selected_names,collapse="+"),")" ))
selected_covs_selection<-list()
selected_covs_outcome<-list()
selected_covs_selection[[1]]<-c("treat:PERS_INC3","treat:MOS_AFDC8")
selected_covs_outcome[[1]]<-c("FEMALE","BLACK")
selected_covs_selection[[2]]<-c("treat:EARN_YR","treat:R_HOME1")
selected_covs_outcome[[2]]<-c("AGE","FEMALE","BLACK")
selected_covs_selection[[3]]<-c("treat:EARN_YR","treat:R_HOME1")
selected_covs_outcome[[3]]<-c("AGE","FEMALE","BLACK")

selected_covs_selection[[4]]<-c("treat:EARN_YR","treat:R_HOME1",
                                "treat:WELF5AC","treat:TYPEWORR5")
selected_covs_outcome[[4]]<-c("FEMALE","BLACK","R_HOME1","HH_INC5")

selected_covs_selection[[5]]<-c( "treat:WELF5AC","treat:TYPEWORR5","treat:IMPRCNTR4","treat:DRG_SUMP2","treat:TYPEJOBB9")
selected_covs_outcome[[5]]<-c("FEMALE","BLACK","IMPRCNTR4")

selected_covs_selection[[6]]<-c("treat:IMPRCNTR4","treat:DRG_SUMP2","treat:TYPEJOBB9")
selected_covs_outcome[[6]]<-c("FEMALE","BLACK")

library(doParallel)
library(foreach)
cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))




for (i in c(3:6)) {
  # prepare data
  week<-selected_weeks[i]
  print (paste0("Results for week ", week))
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome=logwage_week))
  
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome = logwage_week, Lee_data_all_covariates,
                     Lee_data[,baseline_varnames])
  leedata_cov[is.na( leedata_cov)]<-0
  
  leedata_week$weights<-Lee_data$DSGN_WGT.y
  leedata_cov$weights<-Lee_data$DSGN_WGT.y
  min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
  max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  
  ## monotonicity-preserving bounds ##
  glm.fit<-estimate_selection(form=form_nonmonotone_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                              names_to_include = selected_covs_selection[[i]] )
  s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  mean(p.0.star<1)


  
 # estimates_nonmonotone[,i]<-GetBounds(leebounds_wout_monotonicity(leedata_cov,p.0.star))
  
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat= s.hat,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=selected_covs_outcome[[i]],min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,ortho=TRUE,c_quant=0,weight=Lee_data$DSGN_WGT.y) 
 
   
  orthoestimates_postlasso[,i]<-GetBounds(leebounds_ortho_result)


  estimated_orthobounds_bb<-main_bb(leedata_cov,N_rep=N_rep,function_name=second_stage_wrapper,
                                    y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat,
                                    inds_helps=leebounds_ortho_result$inds_helps,weight=Lee_data$DSGN_WGT.y)
  CR_ortho[,i]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb,ATE_est=    orthoestimates_postlasso[,i],ci_alpha=ci_alpha)
  IM_ortho[,i]<-imbens_manski(estimated_orthobounds_bb,orthoestimates_postlasso[,i], ci_alpha=ci_alpha)
  if (FALSE) {
  res<-summary_subjects_positive_lower_bound(leedata_cov_total=leedata_cov[,c("treat","selection","outcome",baseline_varnames)],
                                             s.hat=s.hat,y.hat=leebounds_ortho_result$y.hat,
                                             
                                             quantile_grid_size=quantile_grid_size,
                                             variables_for_outcome=baseline_varnames,
                                             form_outcome=paste0("outcome~",paste0(baseline_varnames,collapse="+")),weights=Lee_data$DSGN_WGT.y)
  estimates_plb[,i]<-GetBounds(res)
  frac_positive[i]<-GetFraction(res)
  }
}

## positive lower bound
estimates_table<-rbind(orthoestimates_postlasso)
estimates_table<-t(estimates_table)
colnames(estimates_table)<-c("Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone")

CR_table<-rbind(CR_ortho)
CR_table<-t(CR_table)
colnames(CR_table)<-c("Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone")


IM_table<-rbind(IM_ortho)
IM_table<-t(IM_table)
estimates_table<-apply(estimates_table,2,round,3)
CR_table<-apply(CR_table,2,round,3)
IM_table<-apply(IM_table,2,round,3)

print("Saving estimates in STEP3_Estimate_Bounds/csv/ ...")
sink(file=NULL)
closeAllConnections()
write.csv(estimates_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_estimates_lasso.csv"))
write.csv(CR_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_CR_lasso.csv"))
write.csv(IM_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_IM_lasso.csv"))

