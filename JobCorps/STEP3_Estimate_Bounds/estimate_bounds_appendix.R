rm(list=ls())

## load packages
print ("Loading packages ...")
library(sas7bdat)
library(stats)
library(expm)
library(feather)
library(reldist)
library(xtable)
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
selected_weeks<-90
N_rep=300
ci_alpha=0.05
quantile_grid_size=0.01

# Column 1
estimates_nonmonotone<-matrix(0,2,length(selected_weeks))
CR_nonmonotone<-matrix(0,2,length(selected_weeks))
# Column 2
orthoestimates<-matrix(0,2,length(selected_weeks))
CR_ortho<-matrix(0,2,length(selected_weeks))
IM_ortho<-matrix(0,2,length(selected_weeks))
# Column 3
orthoestimates2<-matrix(0,2,length(selected_weeks))
CR_ortho2<-matrix(0,2,length(selected_weeks))
IM_ortho2<-matrix(0,2,length(selected_weeks))


#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

selected_names<-unique(setdiff(c( "DRG_SUMP2","IMP_PRO1","FRQ_POT3","R_HOME1","TYPEWORR5","MOS_AFDC8",colnames(Lee_data_covariates)),c("treat","selection","outcome","(Intercept)","X.Intercept.")))
form_nonmonotone<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))
form_nonmonotone2<-as.formula(paste0("selection~(treat+AGE+EARN_YR+WKEARNR+HRSWK_JR+MOSINJOB)*(", paste0(selected_names,collapse="+"),")"))


i=1
week<-selected_weeks[i]
print (paste0("Results for week ", week))
  
hwh_name<-paste0("HWH",week)
earn_name<-paste0("EARNH",week)
  
logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0
  
leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome=logwage_week))
  
leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome = logwage_week, Lee_data_all_covariates,Lee_data_covariates[,baseline_varnames])
leedata_cov[is.na( leedata_cov)]<-0
  
leedata_week$weights<-Lee_data$DSGN_WGT.y
leedata_cov$weights<-Lee_data$DSGN_WGT.y
min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  
  
  
### COLUMN 1
print ("Computing  Lee (2009) estimates without monotonicity (column 1) ...")
glm.fit<-glm(form=form_nonmonotone,data= leedata_cov[,setdiff(colnames(leedata_cov),"outcome")],family="binomial")
## reporting first stage coefficients

glm.fit$coefficients

mytable<-sapply(glm.fit$coefficients[glm.fit$coefficients!=0],round,3)[c(3:35,1,36:68,2)]
mytable_new<-matrix("0",34,3)
mytable_new[,1]<-names(mytable)[1:34]
#mytable_new[,3]<-names(mytable)[38:74]

mytable_new[,2]<-mytable[1:34]
mytable_new[,3]<-mytable[35:68]

res<-predict_selection(glm.fit, leedata_cov)
s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
inds_helps<-(p.0.star<=1)
inds_hurts<-(!inds_helps)

training_data<-leedata_cov[leedata_cov$selection == 1 & leedata_cov$treat==1 , ]
myweights<-leedata_cov$weights[leedata_cov$selection == 1 & leedata_cov$treat==1]
q_model_helps<-quantreg::rq(outcome~.,data=training_data[,c("outcome","DRG_SUMP2","IMP_PRO1","FRQ_POT3","R_HOME1","TYPEWORR5","MOS_AFDC8",baseline_varnames)],tau=0.97,weights=myweights )

training_data<-leedata_cov[leedata_cov$selection == 1 & leedata_cov$treat==0 , ]
myweights<-leedata_cov$weights[leedata_cov$selection == 1 & leedata_cov$treat==0 ]
q_model_hurts<-quantreg::rq(outcome~.,data=training_data[,c("outcome","DRG_SUMP2","IMP_PRO1","FRQ_POT3","R_HOME1","TYPEWORR5","MOS_AFDC8",baseline_varnames)],tau=0.95,weights=myweights )


mytable_qr_helps<-sapply(q_model_helps$coefficients,round,3)
mytable_qr_hurts<-sapply(q_model_hurts$coefficients,round,3)

mytable_qr_new<-matrix("0",nrow=34,ncol=3)
mytable_qr_new[,1]<-names(q_model_helps$coefficients)
mytable_qr_new[,2]<-c(mytable_qr_hurts[2:34],mytable_qr_hurts[1])
mytable_qr_new[,3]<-c(mytable_qr_helps[2:34],mytable_qr_helps[1])

mytable_qr_new<-cbind(mytable_new, mytable_qr_new[,2:3])
print(xtable(mytable_qr_new,align="cc|c|c|c|c|" ), type="latex",include.rownames =TRUE)
