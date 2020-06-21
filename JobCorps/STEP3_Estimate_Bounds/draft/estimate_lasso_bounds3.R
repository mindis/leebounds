

i=6
selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula("selection~treat+TSTAYCAT1+PLACEIPC5+HOWSPOKE2+TALK_CW1+R_GETGED1+R_COMM1+R_HOME1+PERS_INC1+PERS_INC2+HH_INC4+HH_INC5+MOS_ANYW12+MOS_FS12+MOS_AFDC12+HRWAGER6+OCC_R2+OCC_R9+REC_JOB2+EVWORKB1+NUMBJOBS3+REASED_R7+HS_D1+OCC_MOTH5+HH144+AGE16+AGE17+JCMSA2+WORK252+WORK3023+WORK344+WORK461+WORK5212+HGC+HGC_FATH+FEMALE+BLACK+CURRJOB.y+MOSINJOB+YR_WORK+treat:TOTALHRS1+treat:TALK_TCH1+treat:ENCR_PAR1+treat:E_READ1+treat:R_HOME1+treat:HEAR_JC5+treat:FRQ_COKE3+treat:FRQ_CIG4+treat:HEALTH2+treat:MARRIAGE1+treat:HH144+treat:NTV_LANG2+treat:NTV_LANG3+treat:SCHL11+treat:WORK313.1+treat:WORK414.1+treat:NHRSED_R+treat:OTHERRAC+TOTALHRS1+TYPEJOBB11+E_MATH1+R_HOME1+FRQ_POT2+FRQ_ALC4+PY_POT1+EV_ALCHL1+HEALTH3+HH_INC1+GOTFS11+GOTAFDC11+REASLEFT3+NCHLD3+WELF_KID2+WELF_KID4+HH141+WORK84+WORK281B+WORK292+WORK482+WELF121B+treat:TYPEWORR5+treat:CPAROLE21+treat:MARRCAT11+treat:WELF_KID2+treat:SCHL253+treat:WORK2C+TRADWANT1+TALKSTAY1+TALK_ADL1+E_ALONG1+KNEW_JC1+DRG_SUMP1+PERS_INC4+TYPEED_R6+HOUS_ARR2+OCC_FATH8+F_WORK141+WELF_KID1+WORK83+WORK422+WORK501+WORK511+WELF9AB3+EARN_CMP+treat:N_GUILTY6+treat:HRWAGER5+treat:MOSINJOB5+treat:WORK41+treat:WORK422+treat:WELF5AC+treat:WELF12ABC+IMP_CW1+IMPRCNTR4+GOVPRG_R1+HOUS_ARR1+SCHL1312+SCHL153+WORK26+WELF3123+NUMB_HH+AGE+treat:PLACEIPC6+treat:IMP_PRO1+treat:TALK_PRO1+treat:TALK_FRD1+treat:HEALTH4+treat:EVWORKB1+treat:VOC_D1+treat:HS_D1+treat:HH145+treat:JCMSA2+treat:WORK301+TOTALHRS4+IMP_PRO1+TYPEJOBB7+TYPEJOBB9+MOSTIMPR1+MOSTIMPR6+OTHMISC1+EVARRST11+MOS_TRTR1+PY_CIG1+YR_WORK11+REASED_R4+HHMEMB5+OCC_FATH7+M_WORK141+SCHL243+WORK93+WORK182+WORK222+WORK421+WORK521+WELF10AC+WELF121BC+treat:TALKSTAY1+treat:TYPEJOBB9+treat:MOSTIMPR1+treat:INFO_JC3+treat:FRQ_POT3+treat:FRQ_ALC4+treat:HH_INC5+treat:MOS_FS8+treat:REC_JOB1+treat:YR_WORK11+treat:REASLEFT3+treat:M_WORK141+treat:HH141+treat:SCHL202+treat:WORK462+treat:MOSINJOB+PLACEIPC6+ENCR_CW1+E_SPCJOB1+EV_CIG1+HEALTH4+REASLEFT1+WELF_KID3+JCMSA3+SCHL233+SCHL401B+SCHL5012+TRNG322+WORK54+WORK202+WORK37123+WORK481+WORK491+HGC_MOTH+treat:TELEMODE1+treat:IMPRCNTR4+treat:DRG_SUMP1+treat:PERS_INC2+treat:OCC_R5+treat:MOSINJOB8+treat:SCHL4912+treat:WORK3912+treat:WELF412")

selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))


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

## estimate selection
glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std, leedata=leedata_cov,selection_function_name="rlassologit",
                            penalty=list(lambda=300))
res<-predict_selection(glm.fit,leedata_cov)
s.hat = data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)

p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
inds_helps<-(p.0.star<=1)
inds_hurts<-(!inds_helps)
mean(p.0.star<1)

## inds_helps
cov_names<-c("DRG_SUMP2","IMPRCNTR4","TYPEJOBB9")
group_weight<-1:length(cov_names)
mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
mygroup<-as.numeric(mygroup)

mygroup<-group_by(data.frame(group=mygroup),group) %>%
  count %>%
  inner_join(data.frame(group=mygroup,MPRID=Lee_data$MPRID,Lee_data_all_covariates[,cov_names]))
## unite groups of size less than 30 into 1 group
mygroup$group[mygroup$n<=20]<-1
mygroup<-ungroup(mygroup)
mygroup$group[mygroup$group %in% c(0,1,2,3)]<-0
### what are important covs?

