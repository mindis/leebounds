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
#sink(paste0(my_path,"/JobCorps/STEP3_Estimate_Bounds/Table1_Col123.log"))
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
N_rep=300
ci_alpha=0.05
quantile_grid_size=0.01

orthoestimates_postlasso<-matrix(0,2,length(selected_weeks))
CR_ortho<-matrix(0,2,length(selected_weeks))
IM_ortho<-matrix(0,2,length(selected_weeks))

estimates_plb<-matrix(0,2,length(selected_weeks))
frac_positive<-rep(0,length(selected_weeks))
#### Estimate selection equation
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]

selected_covs_selection<-list()
selected_covs_outcome<-list()

penalty<-rep(0,6)
penalty[2]<-850
penalty[3]<-550
penalty[6]<-550
my_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept.","MPRID" ))
form_nonmonotone_lasso_std<-as.formula(paste0("selection~(treat)*(", paste0(my_names,collapse="+"),")"))
covs<-list()
prob_helps<-rep(0,6)

selected_covs_selection<-list()
selected_covs_selection[[1]]<-c("treat:EARN_YR","treat:R_HOME1")
selected_covs_selection[[2]]<-c("treat:EARN_YR","treat:R_HOME1")
selected_covs_selection[[3]]<-c("treat:EARN_YR","treat:R_HOME1","treat:AGE")
selected_covs_selection[[4]]<-c("treat:EARN_YR","treat:R_HOME1","treat:AGE","treat:FEMALE")
selected_covs_selection[[5]]<-c(baseline_varnames,paste0("treat:",baseline_varnames),"treat:EARN_CMP")
selected_covs_selection[[6]]<-unique(c(baseline_varnames,"treat:EARN_YR","treat:EARN_CMP"  ))
selected_covs_outcome<-list()

frac_plb<-rep(0,6)
orthoestimates_plb<-matrix(0,2,6)
CR_ortho_plb<-matrix(0,2,6)
for (i in 2:6) {
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
  glm.fit<-estimate_selection(form=form_nonmonotone_lasso_std,leedata=leedata_cov,selection_function_name = "rlassologit",
                              penalty=list(lambda=300)
  )
  covs[[i]]<-names(glm.fit$coefficients)
  
  selected_names<-setdiff( c(setdiff(covs[[i]],grep("treat:",covs[[i]],value=TRUE)),unlist(strsplit(grep("treat:",covs[[i]],value=TRUE),"treat:"))),   c("","treat","(Intercept)") )
  
  form_nonmonotone_ss<-as.formula(paste0("selection~treat*(",paste0(selected_names,collapse="+"), ")*(",paste0(selected_names,collapse="+"),")" ))
  
  
  
  glm.fit<-estimate_selection(form=form_nonmonotone_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                              names_to_include=selected_covs_selection[[i]],
                              penalty=list(lambda=850))
  
  s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(!inds_helps)
  prob_helps[i]<-mean(p.0.star<1)
  
  lm.fit<-rlasso(as.formula('outcome~.'), data=leedata_cov[leedata_cov$selection==1 & leedata_cov$treat==1,setdiff(colnames(leedata_cov),
                                                                                                                   c("selection","MPRID","weights" ))])
  
  selected_covs_outcome[[i]]<-setdiff(names(lm.fit$coefficients)[lm.fit$coefficients!=0],c("(Intercept)"))
  # estimates_nonmonotone[,i]<-GetBounds(leebounds_wout_monotonicity(leedata_cov,p.0.star))
  if (i==4) {
    selected_covs_outcome[[i]]<-c(selected_covs_outcome[[i]], "BLACK")
  }
  if (i==5) {
    selected_covs_outcome[[i]]<-unique(c(selected_covs_outcome[[i]], baseline_varnames))
  }
  leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat= s.hat,
                                          quantile_grid_size = quantile_grid_size,
                                          variables_for_outcome=setdiff(selected_covs_outcome[[i]], c("REC_ED5","REC_ED8")),min_wage=min_wage,
                                          max_wage=max_wage,distribution_functionname="rq",
                                          sort_quantiles= TRUE,ortho=TRUE,c_quant=0,weight=Lee_data$DSGN_WGT.y) 
  

  if (i==2) {
    selected_vars<-c("FEMALE","EARN_YR","EARN_CMP","WKEARNR","R_HOME1","NUMBJOBS","MONINED"
    )
  }
  if (i==3) {
    selected_covs_outcome[[i]]<-unique( selected_covs_outcome[[i]],c("PAY_RENT1", "HH_INC5" ,  "PERS_INC1", "HRWAGER"  , "WKEARNR" ,  "FEMALE" ))
  }
  
  if (i==4) {
    selected_vars<-baseline_varnames
    
   
  }
  if (i == 5) {
    selected_vars<-unique(c("FEMALE","EARN_YR","EARN_CMP","WKEARNR","R_HOME1","NUMBJOBS","MONINED","BLACK"))
    
  }
     
  if (i == 6) {
    selected_vars<-c("FEMALE","EARN_YR","EARN_CMP","WKEARNR","NUMBJOBS","MONINED")
    
  }  

    is_positive<-summary_subjects_positive_lower_bound2(leedata_cov_total=leedata_cov,
                                             s.hat=   leebounds_ortho_result$s.hat,y.hat=  leebounds_ortho_result$y.hat,
                                             weight=Lee_data$DSGN_WGT.y,myfun=lm,form_outcome=as.formula(paste0("outcome~",paste0(selected_vars,collapse="+")) ))
  
    frac_plb[i]<-mean(is_positive)
 leebounds_ortho_result2<-second_stage_wrapper(leedata=leedata_cov[  is_positive,],
                                                inds_helps=inds_helps[is_positive],
                                                s.hat= leebounds_ortho_result$s.hat[  is_positive,],
                                                y.hat=leebounds_ortho_result$y.hat[is_positive,],
                                                weight=Lee_data$DSGN_WGT.y[  is_positive],ortho=TRUE) 
 orthoestimates_plb[,i]<-GetBounds(leebounds_ortho_result2)
 
 if (i==3) {
   leebounds_ortho_result2<-summary_subjects_positive_lower_bound(leedata=leedata_cov,
                                                 s.hat= leebounds_ortho_result$s.hat,
                                                 y.hat=leebounds_ortho_result$y.hat,form_outcome=as.formula(paste0("outcome~",paste0(selected_vars,collapse="+"))),
                                                 weight=Lee_data$DSGN_WGT.y) 
   orthoestimates_plb[,i]<-GetBounds(leebounds_ortho_result2)
 }
 estimated_orthobounds_bb<-main_bb(leedata_cov,N_rep=N_rep,function_name=summary_subjects_positive_lower_bound,
                                   y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat
                                    ,weight=Lee_data$DSGN_WGT.y,form_outcome=as.formula(paste0("outcome~",paste0(selected_vars,collapse="+"))))
 CR_ortho[,i]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb[!is.na(apply( estimated_orthobounds_bb,1,sum)),],ATE_est=    orthoestimates_plb[,i],ci_alpha=ci_alpha)
 IM_ortho[,i]<-imbens_manski(estimated_orthobounds_bb,orthoestimates_plb[,i], ci_alpha=ci_alpha)
 
 
}
