rm(list=ls())

## load packages
print ("Loading packages ...")
library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)
library(hdm)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
#sink(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/Table1_Col4_lasso.log"))
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
selected_weeks<-c(45,90,104,135,180,208)
weeks<-1:208
N_rep=300
ci_alpha=0.05
quantile_grid_size=0.01

orthoestimates_plb<-matrix(0,2,length(weeks))
frac_positive<-rep(0,length(weeks))
#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

selected_covs_selection<-list()
selected_covs_outcome<-list()

my_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept.","MPRID" ))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(my_names,collapse="+"),")"))
covs<-list()
prob_helps<-rep(0,6)


is_positive_list<-list()
treatment_effect_lower_bound_list<-list()

### step 1: find subjects with positive lower bound in each week
for (i in 80:120) {
  week<-weeks[i]
  
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
  glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function_name = "rlassologit",
                              penalty=list(lambda=300)
  )
  covs[[i]]<-names(glm.fit$coefficients)
  
  selected_names<-setdiff( c(setdiff(covs[[i]],grep("treat:",covs[[i]],value=TRUE)),unlist(strsplit(grep("treat:",covs[[i]],value=TRUE),"treat:"))),   c("","treat","(Intercept)") )
  
  form_nonmonotone_ss<-as.formula(paste0("selection~treat*(",paste0(selected_names,collapse="+"), ")*(",paste0(selected_names,collapse="+"),")" ))
  
  
  glm.fit<-estimate_selection(form=form_nonmonotone_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                              names_to_include=c("treat:EARN_YR","treat:R_HOME1"),
                              penalty=list(lambda=850))
  
  s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  prob_helps[i]<-mean(p.0.star<1)
  
  leebounds_ortho_result<-first_stage_wrapper(leedata_cov,weights=Lee_data$DSGN_WGT.y,
                                              variables_for_outcome=c("EARN_YR","EARN_CMP","FEMALE","BLACK","PERS_INC1","PERS_INC4","HH_INC5","WKEARNR"),
                                              quantile_grid_size=0.01,
                                              sort_quantiles=TRUE,
                                              s.hat=s.hat)
  
  selected_vars<-c(baseline_varnames,"EARN_CMP")
  res<-treatment_effect_lower_bound(leedata_cov_total=leedata_cov[,c("treat","selection","outcome",selected_vars)],
                                             s.hat=leebounds_ortho_result$s.hat,y.hat=leebounds_ortho_result$y.hat,
                                             
                                             quantile_grid_size=quantile_grid_size,
                                             variables_for_outcome=baseline_varnames,
                                             form_outcome=paste0("outcome~",paste0(selected_vars,collapse="+")),weights=Lee_data$DSGN_WGT.y,
                                             cutoff_positive=0.001)
  
  treatment_effect_lower_bound_list[[i]]<-res$treat_effect_lower_bound
}
## define target population of interest: positive lower bound at any week in (45,90,104,135,180,208)
ps_any_week<-apply(matrix(unlist(treatment_effect_lower_bound_list),nrow=9145),1,mean)>0
for (i in 5:208) {
  
  week<-weeks[i]
  
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
  glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function_name = "rlassologit",
                              penalty=list(lambda=300)
  )
  covs[[i]]<-names(glm.fit$coefficients)
  
  selected_names<-setdiff( c(setdiff(covs[[i]],grep("treat:",covs[[i]],value=TRUE)),unlist(strsplit(grep("treat:",covs[[i]],value=TRUE),"treat:"))),   c("","treat","(Intercept)") )
  
  form_nonmonotone_ss<-as.formula(paste0("selection~treat*(",paste0(selected_names,collapse="+"), ")*(",paste0(selected_names,collapse="+"),")" ))
  
  
  #selected_covs_selection[[i]]<-unique(c(baseline_varnames,paste0("treat:",baseline_varnames),"treat:EARN_CMP","treat:R_HOME1" ))
  glm.fit<-estimate_selection(form=form_nonmonotone_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                              names_to_include=unique(c(baseline_varnames,paste0("treat:",baseline_varnames) )),
                              penalty=list(lambda=850))
  
  s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  prob_helps[i]<-mean(p.0.star<1)
  
  ### first stage is estimated on the full sample
  fs_result<-first_stage_wrapper(leedata_cov,weights=Lee_data$DSGN_WGT.y,
                                              variables_for_outcome=c("EARN_YR","EARN_CMP","FEMALE","BLACK","PERS_INC1","PERS_INC4","HH_INC5","WKEARNR"),
                                              quantile_grid_size=0.01,
                                              sort_quantiles=TRUE,
                                              s.hat=s.hat)
  ## bounds (moment functions) are evaluated for a subset ps_any_week
  leebounds_ortho_result<-second_stage_wrapper(leedata=leedata_cov[  ps_any_week ,],
                                                inds_helps=inds_helps[ps_any_week],
                                                s.hat=  fs_result$s.hat[  ps_any_week,],
                                                y.hat= fs_result$y.hat[ps_any_week,],
                                                weight=Lee_data$DSGN_WGT.y[ps_any_week],ortho=TRUE,c_quant=0) 
  orthoestimates_plb[,i]<-GetBounds(leebounds_ortho_result)
  ## week-specific bounds
  
  
}
selected_weeks<-5:208
Figure3_dataset<-data.frame(weeks=weeks[selected_weeks],
                            lb_positive=orthoestimates_plb[1,selected_weeks],
                            ub_positive=orthoestimates_plb[2,selected_weeks]
                            
)
write.csv(Figure3_dataset,paste0(my_path,"/JobCorps/STEP5_Print_Figures/csv/Figure3_selected_covs3.csv"))
