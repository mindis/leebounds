## load packages
print ("Loading packages ...")
library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)
library(hdm)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
print ("Loading data ...")
Lee_data_covariates<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps_data/dataLee2009covariates3.feather"))

Lee_data<-as.data.frame(Lee_data_covariates)

source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/orthogonal_correction.R"))
source(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/utils.R"))

N_rep=300
ci_alpha=0.05
quantile_grid_size=0.01


#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

selected_covs_selection<-c("treat:EARN_YR","treat:R_HOME1")


my_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept.","MPRID" ))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(my_names,collapse="+"),")"))
covs<-list()
prob_helps<-rep(0,6)
selected_weeks<-c(45,90,104,135,180,208)
N_rep=300
ci_alpha=0.05
quantile_grid_size=0.01
week<-90
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
i=2

## monotonicity-preserving bounds ##
glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function_name = "rlassologit",
                            penalty=list(lambda=300)
)
covs[[i]]<-names(glm.fit$coefficients)

selected_names<-setdiff( c(setdiff(covs[[i]],grep("treat:",covs[[i]],value=TRUE)),unlist(strsplit(grep("treat:",covs[[i]],value=TRUE),"treat:"))),   c("","treat","(Intercept)") )

form_nonmonotone_ss<-as.formula(paste0("selection~treat*(",paste0(selected_names,collapse="+"), ")*(",paste0(selected_names,collapse="+"),")" ))



glm.fit<-estimate_selection(form=form_nonmonotone_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                            names_to_include=selected_covs_selection,
                            penalty=list(lambda=850))

s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
inds_helps<-(p.0.star<=1)
inds_hurts<-(!inds_helps)
prob_helps[i]<-mean(p.0.star<1)

mytable<-sapply(glm.fit$coefficients[glm.fit$coefficients!=0],round,3)
mytable_new<-matrix("0",11,3)
mytable_new[,1]<-names(mytable)[1:11]
mytable_new[,2]<-mytable[1:11]


lm.fit<-rlasso(as.formula('outcome~.'), data=leedata_cov[leedata_cov$selection==1 & leedata_cov$treat==1,setdiff(colnames(leedata_cov),
                                                                                                                 c("selection","MPRID","weights" ))])

selected_covs_outcome<-setdiff(names(lm.fit$coefficients)[lm.fit$coefficients!=0],c("(Intercept)"))


training_data<-leedata_cov[leedata_cov$selection == 1 & leedata_cov$treat==1 , ]
myweights<-leedata_cov$weights[leedata_cov$selection == 1 & leedata_cov$treat==1]
q_model_helps<-quantreg::rq(outcome~.,data=training_data[,c("outcome",selected_covs_outcome)],tau=0.97,weights=myweights )

training_data<-leedata_cov[leedata_cov$selection == 1 & leedata_cov$treat==0 , ]
myweights<-leedata_cov$weights[leedata_cov$selection == 1 & leedata_cov$treat==0 ]
q_model_hurts<-quantreg::rq(outcome~.,data=training_data[,c("outcome",selected_covs_outcome)],tau=0.95,weights=myweights )

mytable_qr_helps<-sapply(q_model_helps$coefficients,round,3)
mytable_qr_hurts<-sapply(q_model_hurts$coefficients,round,3)

mytable_qr_new<-matrix("0",nrow=7,ncol=3)
mytable_qr_new[,1]<-c(names(q_model_helps$coefficients)[2:7], names(q_model_helps$coefficients)[1])
mytable_qr_new[,2]<-c(mytable_qr_hurts[2:7],mytable_qr_hurts[1])
mytable_qr_new[,3]<-c(mytable_qr_helps[2:7],mytable_qr_helps[1])

mytable_qr_new<-cbind(mytable_new, mytable_qr_new[,2:3])
print(xtable(mytable_qr_new,align="cc|c|c|c|c|" ), type="latex",include.rownames =TRUE)
