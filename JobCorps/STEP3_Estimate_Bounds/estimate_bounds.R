rm(list=ls())

## load packages
print ("Loading packages ...")
library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)
library(reldist)
library(xtable)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
sink(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/Table1_Col123.log"))
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
N_rep=500
ci_alpha=0.05
quantile_grid_size=0.01

# Column 1
estimates_nonmonotone<-matrix(0,2,length(selected_weeks))
CR_nonmonotone<-matrix(0,2,length(selected_weeks))
IM_nonmonotone<-matrix(0,2,length(selected_weeks))
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

selected_names<-unique(setdiff(c( "DRG_SUMP2","TYPEJOBB9","REASED_R4","IMP_PRO1","FRQ_POT3","R_HOME1","WELF5AC","TYPEWORR5","MOS_AFDC8",colnames(Lee_data_covariates)),c("treat","selection","outcome","(Intercept)","X.Intercept.")))
form_nonmonotone<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))
form_nonmonotone2<-as.formula(paste0("selection~(treat+AGE+EARN_YR+WKEARNR+HRSWK_JR+MOSINJOB)*(", paste0(selected_names,collapse="+"),")"))


for (i in c(1:6)) {
  # prepare data
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
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  estimates_nonmonotone[,i]<-GetBounds(leebounds_wout_monotonicity(leedata_cov,p.0.star))
  
  if (TRUE) {
  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,p.0.star=p.0.star,mydata=leedata_week[,c("treat","selection","outcome")],N_rep=N_rep)
  ## confidence region for identified set
  print ("Computing confidence region for Lee (2009) without monotonicity  estimates by regular bootstrap (column 3) ...")
  CR_nonmonotone[,i]<-compute_confidence_region(bounds_bb,estimates_nonmonotone[,i], ci_alpha=ci_alpha )
  IM_nonmonotone[,i]<-imbens_manski(bounds_bb,estimates_nonmonotone[,i], ci_alpha=ci_alpha)
  }
  ### COLUMN 2
  print ("Computing  Lee (2009) estimates without monotonicity (column 2) ...")
  
  
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=c(baseline_varnames),min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,weights=Lee_data$DSGN_WGT.y)
  
  
  orthoestimates[,i]<-GetBounds(leebounds_ortho_result)
  
  
  
  
  if (TRUE) {
  estimated_orthobounds_bb<-main_bb(leedata_cov,N_rep=N_rep,function_name=second_stage_wrapper,
                                      y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat,
                                       inds_helps=leebounds_ortho_result$inds_helps)
  CR_ortho[,i]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb,ATE_est=  orthoestimates[,i],ci_alpha=ci_alpha)
  IM_ortho[,i]<-imbens_manski(estimated_orthobounds_bb,orthoestimates[,i], ci_alpha=ci_alpha)
  
  }
  ### COLUMN 3
  print ("Computing  Lee (2009) estimates without monotonicity (column 3) ...")
  
  
  leebounds_ortho_result2<-ortho_leebounds(leedata_cov=leedata_cov,form=form_nonmonotone2,
                                          selection_function_name="glm",
                                          selection_function=glm,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=c(baseline_varnames),min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,weights=Lee_data$DSGN_WGT.y)
  
  
  orthoestimates2[,i]<-GetBounds(leebounds_ortho_result2)
  
  if (TRUE) {
  estimated_orthobounds_bb2<-main_bb(leedata_cov,N_rep=N_rep,function_name=second_stage_wrapper,
                                    y.hat= leebounds_ortho_result2$y.hat,s.hat=leebounds_ortho_result2$s.hat,
                                    inds_helps=leebounds_ortho_result2$inds_helps)
  CR_ortho2[,i]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb2,ATE_est=  orthoestimates2[,i],ci_alpha=ci_alpha)
  IM_ortho2[,i]<-imbens_manski(estimated_orthobounds_bb2,orthoestimates2[,i], ci_alpha=ci_alpha)
  
  }
}
  


estimates_table<-rbind(estimates_nonmonotone,orthoestimates,orthoestimates2)
estimates_table<-t(estimates_table)
colnames(estimates_table)<-c("Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone", "OrthoLee_28_lb","OrthoLee_28_ub",
                             "OrthoLee_28_lb2","OrthoLee_28_ub2"
                             )
CR_table<-rbind(CR_nonmonotone,CR_ortho,CR_ortho2)
CR_table<-t(CR_table)
colnames(CR_table)<-c("Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone",
                      "OrthoLee_28_lb","OrthoLee_28_ub",
                      "OrthoLee_28_lb2","OrthoLee_28_ub2")


IM_table<-rbind(IM_nonmonotone,IM_ortho,IM_ortho2)
IM_table<-t(IM_table)
estimates_table<-apply(estimates_table,2,round,3)
CR_table<-apply(CR_table,2,round,3)
IM_table<-apply(IM_table,2,round,3)

print("Saving estimates in STEP3_Estimate_Bounds/csv/ ...")
sink(file=NULL)
closeAllConnections()
write.csv(estimates_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col123_estimates.csv"))
write.csv(CR_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col123_CR.csv"))
write.csv(IM_table,paste0(my_path,"JobCorps/STEP5_Print_Tables/csv/Table1_Col123_IM.csv"))

### save as latex table
table<-print_table(estimates_table,CR_table,IM_table,digs=3)
table<-print(xtable(table,type="latex",include.rownames =FALSE ))
write.table(table,paste0(my_path,"/JobCorps/STEP5_Print_Tables/Table123.txt"),append=TRUE)

