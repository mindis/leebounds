rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")
### load libraries
library(expm)
library(feather)
library(hdm)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
sink(paste0(my_path,"/JobCorps/STEP4_MISC/Figure4.log"))
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

source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/orthogonal_correction.R"))




Nboot=300
ci_alpha=0.05
quantile_grid_size=0.01
week1=104
week2=208
## parametrization of half-2d shere
phis<-seq(0,(pi),by=0.1)
### selection equation specification
selected_names<-unique(c("R_GETGED1","R_COMM1","R_HOME1","PERS_INC1","PERS_INC2","HH_INC4",
                  "HH_INC5","MOS_ANYW12"," MOS_FS12","MOS_AFDC12","EVWORKB1",
                  "REASED_R7","HS_D1","HH144",   "AGE16",     "HGC"    ,    "FEMALE" ,    "BLACK" ,    
                   "MOSINJOB",   "YR_WORK",    "TYPEJOBB11", "E_MATH1" ,  "FRQ_POT2" ,  "FRQ_ALC4",   "EV_ALCHL1" , "HH_INC1" ,
                  "REASLEFT3", "WELF_KID2",  "WELF_KID4" , "WORK292" ,   "TRADWANT1" , "TALKSTAY1",  "TALK_ADL1",
                  "E_ALONG1" ,  "KNEW_JC1",   "DRG_SUMP1" ,
                  "PERS_INC4",  "TYPEED_R6",  "HOUS_ARR2" , "OCC_FATH8",  "F_WORK141" , "WELF_KID1" , "WORK83",
                  "WORK422" ,   "WORK501"  , 
                   "WORK511" ,   "WELF9AB3",   "EARN_CMP" ,  "GOVPRG_R1" , "WORK2C"  ,   "N_GUILTY6" , "HRWAGER5",
                  "MOSINJOB5",  "WORK41"   , 
                 "WELF5AC" ,   "WELF12ABC" ,baseline_varnames
                  ))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))
#form_nonmonotone<-as.formula("selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412")


############### ESTIMATE SELECTION EQUATION FOR WEEK 104 and 208 ####

hwh_name<-paste0("HWH",week1)
earn_name<-paste0("EARNH",week1)

logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0

selection1<-logwage_week>0
outcome1<-logwage_week

hwh_name<-paste0("HWH",week2)
earn_name<-paste0("EARNH",week2)

logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0

selection2<-logwage_week>0
outcome2<-logwage_week



selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula(paste0("selection~treat*(", paste0(selected_names,collapse="+"),")"))



## leedata_cov
print (paste0("Estimating selection equation for week ",  week1))
leedata_cov1<-cbind(treat=Lee_data$TREATMNT.y,selection=selection1,outcome = outcome1, Lee_data_all_covariates)
leedata_cov1[is.na( leedata_cov1)]<-0

### week 104: do not include baseline_varnames
glm.fit1<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov1,selection_function = glm,selection_function_name = "rlassologit")
res<-predict_selection(glm.fit1, leedata_cov1)
s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
p.0.star1=s.hat$s.0.hat/s.hat$s.1.hat


## leedata_cov
print (paste0("Estimating selection equation for week ",  week2))
leedata_cov2<-cbind(treat=Lee_data$TREATMNT.y,selection=selection2,outcome = outcome2, Lee_data_all_covariates)
leedata_cov2[is.na( leedata_cov2)]<-0

glm.fit2<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov2,selection_function = glm,selection_function_name = "rlassologit")

res<-predict_selection(glm.fit2, leedata_cov2)
s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
p.0.star2=s.hat$s.0.hat/s.hat$s.1.hat

## 
print (paste0("Preparing the set AE^{2+} "))
inds_helps_helps<-(p.0.star1<1)&(p.0.star2<1)
treat_helps=TRUE
form_monotone<-as.formula(paste0("selection~treat+(", paste0(selected_names,collapse="+"),")"))

estimated_orthobounds_week104208<-matrix(0,nrow=2,length(phis))
estimated_orthobounds_CI_week104208<-matrix(0,nrow=2,length(phis))
estimated_orthobounds_bb<-array(0,c(Nboot,2,length(phis)))

leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=(selection1==1)*(selection2==1),Lee_data_all_covariates)
leedata_cov[is.na(leedata_cov)]<-0
leedata_cov<-leedata_cov[inds_helps_helps,]

print (paste0("Estimating selection equqation for subjects in  AE^{2+} "))
glm.fit<-estimate_selection(form=form_monotone,leedata= as.data.frame(leedata_cov[,setdiff(colnames(leedata_cov),"outcome")]),
                           selection_function = rlassologit,
                            selection_function_name = "rlassologit")
res<-predict_selection(glm.fit, leedata_cov)
s.hat = data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)

for (k in 1:length(phis)) {
  if (floor(k/10)==ceiling(k/10)) {
    print (paste0("Estimating support function for angle: 2*pi",round(k/length(phis),2)))
  }

  phi<-phis[k]
  q1<-cos(phi)
  q2<-sin(phi)
  
  
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=(selection1==1)*(selection2==1),outcome=q1*outcome1+q2*outcome2, Lee_data_all_covariates)
  leedata_cov[is.na( leedata_cov)]<-0
  
  leedata_cov<-leedata_cov[inds_helps_helps,]
  
   
  res<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=selected_names,weights=Lee_data$DSGN_WGT.y[inds_helps_helps])
  
  
  estimated_orthobounds_week104208[,k]<-GetBounds(res)
  estimated_orthobounds_bb[,,k]<-main_bb(leedata_cov,N_rep=Nboot,function_name=second_stage_wrapper,
                                               y.hat=res$y.hat,s.hat=res$s.hat,inds_helps=res$inds_helps,weights=Lee_data$DSGN_WGT.y[inds_helps_helps])
  estimated_orthobounds_CI_week104208[,k]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb[,,k],ATE_est=estimated_orthobounds_week104208[,k],ci_alpha=0.05)
  # may need a different critical value if we want uniform CI
}



print (paste0("Best Circle Approximation to Identified Set "))
## R is close to 0.11858 times (1,1,1)
## Identified set is a circle with center beta0 = (beta10,beta20)

Rvec=(estimated_orthobounds_week104208[2,]-estimated_orthobounds_week104208[1,])/2
##
print ("Estimated Radius of Circle Approximation to Identified Set ")
R=mean(Rvec)
print(R)
print ("Estimated Center of Circle Approximation to Identified Set ")
beta10=estimated_orthobounds_week104208[2,1]-R
beta20=estimated_orthobounds_week104208[2,length(phis)/2]-R
print(c(round(beta10,3),round(beta20,3)))
print ("Estimated Radius of Approximation to Confidence Region ")


### Pointwise Confidence Interval for 
Rvec_pointwise=(estimated_orthobounds_CI_week104208[2,]-estimated_orthobounds_CI_week104208[1,])/2
R_pointwise=mean(Rvec_pointwise)
print(R_pointwise)
print ("Estimated Center of Circle Approximation to Confidence Region ")
beta10_pointwise=estimated_orthobounds_CI_week104208[2,1]-R_pointwise
beta20_pointwise=estimated_orthobounds_CI_week104208[2,length(phis)/2]-R_pointwise

print(c(round(beta10_pointwise,3),round(beta20_pointwise,3)))
## rho(0) estimate
print (paste0("Estimated untreated log (wage) for week ", week1))
week1_untreated=weighted.mean(leedata_cov1$outcome[leedata_cov2$selection==1 & leedata_cov1$selection==1 & inds_helps_helps],
              w=Lee_data$DSGN_WGT.y[leedata_cov2$selection==1 & leedata_cov1$selection==1 & inds_helps_helps])
week1_untreated
print (paste0("Estimated untreated log (wage) for week ", week2))
week2_untreated=weighted.mean(leedata_cov2$outcome[leedata_cov2$selection==1 & leedata_cov1$selection==1 & inds_helps_helps],
                                   w=Lee_data$DSGN_WGT.y[leedata_cov2$selection==1 & leedata_cov1$selection==1 & inds_helps_helps])
print (paste0("Estimated untreated log (wage) growth from  ", week1," to ", week2))
week2_untreated
week2_untreated-week1_untreated
sink(file=NULL)
closeAllConnections()