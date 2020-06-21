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
source(paste0(my_path,"/R/orthogonal_correction.R"))
weeks=1:208
selected_weeks<-c(45,90,104,135,180,208)
N_rep=1000
ci_alpha=0.05
quantile_grid_size=0.01

selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))

orthoestimates_postlasso<-matrix(0,2,length(selected_weeks))
CR_ortho_postlasso<-matrix(0,2,length(selected_weeks))

selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."
                                                           ))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))

## save the results
s.hat.lasso_std<-list()
y.hat.lasso_std<-list()

form_final_selection<-list()
selected_covs_selection<-list()

selected_covariates_outcome<-list()
selected_covariates_outcome[[1]]<-c("FEMALE","BLACK","EARN_YR")
selected_covariates_outcome[[2]]<-c("FEMALE","BLACK","EARN_YR","MOSINJOB")
selected_covariates_outcome[[3]]<-c("FEMALE","BLACK","MOSINJOB")
selected_covariates_outcome[[4]]<-c("FEMALE","BLACK","AGE", "EARN_YR")
#selected_covariates_outcome[[5]]<-c("FEMALE" , "BLACK",   "EARN_YR",   "MOSINJOB", "HH_INC5")
selected_covariates_outcome[[5]]<-c("FEMALE" , "BLACK",   "EARN_YR","MOSINJOB" )
selected_covariates_outcome[[6]]<-c("FEMALE","BLACK")
for (i in 1:length(selected_weeks)) {
  
  print(i)
  week<-selected_weeks[i]
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome=logwage_week))
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,outcome = logwage_week, Lee_data_all_covariates)
  leedata_cov[is.na(leedata_cov)]<-0
  min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
  max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
  
 if (i<=3) {
   glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function = rlassologit,selection_function_name = "rlassologit",
                    penalty=list(lambda=400))
   
 } else {
   glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function = rlassologit,selection_function_name = "rlassologit",
                             weights=Lee_data$DSGN_WGT.y,
                             penalty=list(lambda=100),
                             names_to_include=c("IMPRCNTR4","SCHL4912","WELF412","WORK3912","OCC_R5","MOSINJOB8","PERS_INC2"))
   #"WORK3912","WELF412","VOC_D1","OCC_R5","PERS_INC2"c("HS_D1","TELEMODE1","DRG_SUMP1","REASLEFT3","MOSINJOB8","IMPRCNTR4",
  # "SCHL4912")  names_to_include=    c("JCMSA2","WORK301","HH145","HS_D1","HEALTH4","EVWORKB1","TALK_PRO1","IMP_PRO1","PLACEIPC6",
  # "MOSINJOB8")  
   names(glm.fit$coefficients)[order(abs(glm.fit$coefficients))]
 }
    
 
  
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  s.hat.lasso_std[[i]]<-s.hat
  p.0.star<-s.hat.lasso_std[[i]]$s.0.hat/s.hat.lasso_std[[i]]$s.1.hat
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat.lasso_std[[i]],
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=selected_covariates_outcome[[i]],min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,weights=Lee_data$DSGN_WGT.y,c_quant=0) 
  
  orthoestimates_postlasso[,i]<-GetBounds( leebounds_ortho_result)
  
  if (FALSE) {
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
}
orthoestimates_postlasso<-apply(orthoestimates_postlasso,2,round,3)
CR_ortho_postlasso<-apply(  CR_ortho_postlasso,2,round,3)
write.csv(orthoestimates_postlasso,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col5_estimates.csv"))
write.csv(  CR_ortho_postlasso,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col5_CR.csv"))
