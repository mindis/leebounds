rm(list=ls())

## load packages
print ("Loading packages ...")
library(sas7bdat)
library(SDMTools)
library(expm)
library(feather)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
### load data
sink(paste0(my_path,"/JobCorps/STEP2_Estimate_Bounds/Table1_Col1234.log"))
print ("Loading data ...")
Lee_data_covariates<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data_covariates)


source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))

selected_weeks<-c(45,90,104,135,180,208)
N_rep=1000
ci_alpha=0.05
quantile_grid_size=0.01
# Column 1
estimates<-matrix(0,2,length(selected_weeks))
CR<-matrix(0,2,length(selected_weeks))

# Column 2
orthoestimates_monotone<-matrix(0,2,length(selected_weeks))
CR_orthomonotone<-matrix(0,2,length(selected_weeks))

# Column 3
estimates_nonmonotone<-matrix(0,2,length(selected_weeks))
CR_nonmonotone<-matrix(0,2,length(selected_weeks))

# Column 4
orthoestimates<-matrix(0,2,length(selected_weeks))
CR_ortho<-matrix(0,2,length(selected_weeks))

#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]
selected_names<-setdiff(colnames(Lee_data_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept.","FEMALE","BLACK","AGE"))
form_nonmonotone<-as.formula(paste0("selection~(treat+FEMALE+BLACK+AGE)*(", paste0(selected_names,collapse="+"),")"))
form_monotone<-as.formula(paste0("selection~treat+", paste0(baseline_varnames,collapse="+")))


for (i in 1:length(selected_weeks)) {
  # prepare data
  week<-selected_weeks[i]
  print (paste0("Results for week ", week))
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  
  leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome = logwage_week, Lee_data_covariates)
  leedata_cov[is.na( leedata_cov)]<-0

  leedata_week$weights<-Lee_data$DSGN_WGT.y
  
  
  ### COLUMN 1
  
  print ("Computing  Lee (2009) estimates (column 1) ...")
  ## compute bounds 
  estimates[,i]<-GetBounds(leebounds_unknown_sign(leedata_week))
  ##  bootstrap draws of bounds
  print ("Computing  Lee (2009) confidence region by regular bootstrap (column 1) ...")
  #bounds_bb<-main_bb(function_name=leebounds_unknown_sign,mydata=leedata_week,N_rep=N_rep)
  ## confidence region for identified set
  #CR[,i]<-compute_confidence_region(bounds_bb,estimates[,i], ci_alpha=ci_alpha )
  
  ### COLUMN 2 Incorrectly assuming monotonicity 
  ### Adding covariates
  leedata_cov$weights<-Lee_data$DSGN_WGT.y
  glm.fit<-glm(form=form_monotone,data= leedata_cov[,setdiff(colnames(leedata_cov),"outcome")],family="binomial",weights=weights)
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  
  if (week>=90) {
    treat_helps<-TRUE
    p.0.star<-sapply(p.0.star,max,0.99)
  } else {
    treat_helps<-FALSE
    p.0.star<-sapply(p.0.star,min,1.001001)
  }
  res<-ortho_leebounds(leedata_cov,    
                       s.hat=s.hat,
                       quantile_grid_size = quantile_grid_size,
                       variables_for_outcome=baseline_varnames,sort_quantiles=TRUE,weights=Lee_data$DSGN_WGT.y)

  orthoestimates_monotone[,i]<-GetBounds(res)
  
#  estimated_orthobounds_bb<-weighted_bb(leedata_cov,B=N_rep,function_name=second_stage_wrapper, inds_helps=res$inds_helps,
#                                        y.hat=res$y.hat, s.hat=res$s.hat)
#  CR_orthomonotone[,i]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est=  orthoestimates_monotone[,i],ci_alpha=ci_alpha)
  
  
  
  ### COLUMN 3
  print ("Computing  Lee (2009) estimates without monotonicity (column 3) ...")
  glm.fit<-glm(form=form_nonmonotone,data= leedata_cov[,setdiff(colnames(leedata_cov),"outcome")],family="binomial",weights=weights)
  res<-predict_selection(glm.fit, leedata_cov)
  s.hat=data.frame(s.0.hat=res$s.0.hat, s.1.hat=res$s.1.hat)
   p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  estimates_nonmonotone[,i]<-GetBounds(leebounds_wout_monotonicity(leedata_week,p.0.star))
  
#  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,p.0.star=p.0.star,mydata=leedata_week,N_rep=N_rep)
  ## confidence region for identified set
 # print ("Computing confidence region for Lee (2009) without monotonicity  estimates by regular bootstrap (column 3) ...")
 # CR_nonmonotone[,i]<-compute_confidence_region(bounds_bb,estimates_nonmonotone[,i], ci_alpha=ci_alpha )
  

  
  ### COLUMN 4
  print ("Computing  orthogonal estimates with 28 covariates (column 4) ...")
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=baseline_varnames,weights=Lee_data$DSGN_WGT.y )
  
  
  orthoestimates[,i]<-GetBounds(leebounds_ortho_result)
  
  print ("Computing  confidence region for orthogonal estimates with 28 covariates (column 3) by weighted bootstrap...")
  #estimated_orthobounds_bb<-weighted_bb(leedata_cov,B=N_rep,function_name=second_stage_wrapper,
    #                                    y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat,
   #                                     inds_helps=leebounds_ortho_result$inds_helps)
  #CR_ortho[,i]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est=  orthoestimates[,i],ci_alpha=ci_alpha)
  
  
}
estimates_table<-rbind(estimates,orthoestimates_monotone, estimates_nonmonotone,orthoestimates)
estimates_table<-t(estimates_table)
colnames(estimates_table)<-c("Lee_2009_lb","Lee_2009_ub","Lee_2009_lb_cov","Lee_2009_ub_cov","Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone",
                             "OrthoLee_28_lb","OrthoLee_28_ub")
CR_table<-rbind(CR,CR_orthomonotone,CR_nonmonotone,CR_ortho)
CR_table<-t(CR_table)
colnames(CR_table)<-c("Lee_2009_lb","Lee_2009_ub","Lee_2009_lb_cov","Lee_2009_ub_cov","Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone",
                      "OrthoLee_28_lb","OrthoLee_28_ub")

estimates_table<-apply(estimates_table,2,round,3)
CR_table<-apply(CR_table,2,round,3)
print("Saving estimates in STEP2_Estimate_Bounds/csv/ ...")
sink(file=NULL)
closeAllConnections()
write.csv(estimates_table,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col1234_estimates.csv"))
write.csv(CR_table,paste0(my_path,"JobCorps/STEP2_Estimate_Bounds/csv/Table1_Col1234_CR.csv"))
