rm(list=ls())
library(hdm)
library(expm)
library(ranger)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/ABBKK_2002/report_results.R"))

## load data
mydata<-read.csv(paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
dim(mydata)
## clean data
mydata<-mydata[mydata$AGE2<=23 & mydata$SEX_IS_NA == 0, ]

#Parameters

ci_alpha=0.05
## number of bootstrap repetition
N_rep=1000
quantile_grid_size=0.01
exogenous_covariates<-c("AGE2","SEX_NAME","MOM_SCH" ,"MOM_AGE","DAD_SCH" ,"DAD_AGE",
                        "DAREA4",   "DAREA5",   "DAREA6",   "DAREA7",
                        "DAREA11" ,  "DAREA15",
                        "DAREA16" , "DAREA17"  ,"DAREA18" , "DAREA19",
                        
                        "STRATA1",  "STRATA2",  "STRATA3" , "STRATA4", "MOM_AGE_IS_NA","DAD_AGE_IS_NA","MOM_SCH_IS_NA","DAD_SCH_IS_NA"
)
form_nonmonotone<-as.formula("selection ~ treat * (AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                             DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                             STRATA3 + STRATA4 +MOM_AGE_IS_NA+DAD_AGE_IS_NA+MOM_SCH_IS_NA+DAD_SCH_IS_NA )")
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])

## Save Results Column 1
## estimates
estimates<-list()
## confidence region
CR<-list()
estimates_boys<-list()
estimates_girls<-list()

CR_boys<-list()
CR_girls<-list()

estimates_nonmonotone<-list()
CR_nonmonotone<-list()

estimates_nonmonotone_girls<-list()
CR_nonmonotone_girls<-list()

estimates_nonmonotone_boys<-list()
CR_nonmonotone_boys<-list()

orthoestimates<-list()
CR_ortho<-list()

orthoestimates_boys<-list()
CR_ortho_boys<-list()

orthoestimates_girls<-list()
CR_ortho_girls<-list()


### positive lower bound results 

frac_positive<-list()
estimates_plb<-list()
CR_plb<-list()


frac_positive_boys<-list()
estimates_plb_boys<-list()
CR_plb_boys<-list()


frac_positive_girls<-list()
estimates_plb_girls<-list()
CR_plb_girls<-list()


### post-lasso results


orthoestimates_postlasso<-list()
CR_ortho_postlasso<-list()

orthoestimates_postlasso_boys<-list()
CR_ortho_postlasso_boys<-list()

orthoestimates_postlasso_girls<-list()
CR_ortho_postlasso_girls<-list()


### random forest results
orthoestimates_rf<-list()
CR_ortho_rf<-list()

orthoestimates_rf_boys<-list()
CR_ortho_rf_boys<-list()

orthoestimates_rf_girls<-list()
CR_ortho_rf_girls<-list()

leedata<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK)
leedata_boys<-leedata[mydata$SEX_NAME==1,]
leedata_girls<-leedata[mydata$SEX_NAME==0 ,]

############### COLUMN 2: BASIC LEE (2008) BOUNDS UNDER MONOTONICITY ####################

for (subject in c("MATH","READING","WRITING")) {
  leedata$outcome<-mydata[,subject]
  
  ## compute estimbounds 
  estimates[[subject]]<-GetBounds(leebounds_unknown_sign(leedata))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_unknown_sign,mydata=leedata,N_rep=N_rep)
  ## confidence region for identified set
  CR[[subject]]<-compute_confidence_region(bounds_bb,estimates[[subject]], ci_alpha=ci_alpha )
}
p0<-mean(leedata$selection[leedata$treat==0])/mean(leedata$selection[leedata$treat==1])
report_results(estimates,CR)


for (subject in c("MATH","READING","WRITING")) {
  leedata_boys$outcome<-mydata[mydata$SEX_NAME==1,subject]
  
  ## compute estimbounds 
  estimates_boys[[subject]]<-GetBounds(leebounds_unknown_sign(leedata_boys))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_unknown_sign,mydata=leedata_boys,N_rep=N_rep)
  ## confidence region for identified set
  CR_boys[[subject]]<-compute_confidence_region(bounds_bb,estimates_boys[[subject]], ci_alpha=ci_alpha )
}
p0_boys<-mean(leedata_boys$selection[leedata_boys$treat==0])/mean(leedata_boys$selection[leedata_boys$treat==1])


for (subject in c("MATH","READING","WRITING")) {
  leedata_girls$outcome<-mydata[mydata$SEX_NAME==0,subject]
  
  ## compute estimbounds 
  estimates_girls[[subject]]<-GetBounds(leebounds_unknown_sign(leedata_girls))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_unknown_sign,mydata=leedata_girls,N_rep=N_rep)
  ## confidence region for identified set
  CR_girls[[subject]]<-compute_confidence_region(bounds_bb,estimates_girls[[subject]], ci_alpha=ci_alpha )
}
p0_girls<-mean(leedata_girls$selection[leedata_girls$treat==0])/mean(leedata_girls$selection[leedata_girls$treat==1])

################ COLUMN 3: LEE (2008) bounds without mononotonicity and covariates ####################

glm.fit<-estimate_selection(form=form_nonmonotone,selection_function = glm,
                            leedata=leedata_cov,variables_for_selection = setdiff(colnames(leedata_cov),"outcome"))
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat




for (subject in c("MATH","READING","WRITING")) {
  leedata$outcome<-mydata[,subject]
  
  ## compute estimbounds 
  estimates_nonmonotone[[subject]]<-GetBounds(leebounds_wout_monotonicity(leedata,p.0.star=p.0.star))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,mydata=leedata,N_rep=N_rep,p.0.star=p.0.star)
  ## confidence region for identified set
  CR_nonmonotone[[subject]]<-compute_confidence_region(bounds_bb,estimates[[subject]], ci_alpha=ci_alpha )
}
p0<-mean(leedata$selection[leedata$treat==0])/mean(leedata$selection[leedata$treat==1])
report_results(estimates_nonmonotone,CR_nonmonotone)

estimates_boys_nm<-list()
estimates_girls_nm<-list()

CR_boys_nm<-list()
CR_girls_nm<-list()

leedata_boys<-leedata[mydata$SEX_NAME==1,]
for (subject in c("MATH","READING","WRITING")) {
  leedata_boys$outcome<-mydata[mydata$SEX_NAME==1,subject]
  
  ## compute estimbounds 
  estimates_boys_nm[[subject]]<-GetBounds(leebounds_wout_monotonicity(leedata_boys,p.0.star = p.0.star[mydata$SEX_NAME==1]))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,mydata=leedata_boys,N_rep=N_rep,p.0.star =  p.0.star[mydata$SEX_NAME==1])
  ## confidence region for identified set
  CR_boys_nm[[subject]]<-compute_confidence_region(bounds_bb,estimates_boys_nm[[subject]], ci_alpha=ci_alpha )
}

leedata_girls<-leedata[mydata$SEX_NAME==0 ,]
for (subject in c("MATH","READING","WRITING")) {
  leedata_girls$outcome<-mydata[mydata$SEX_NAME==0,subject]
  
  ## compute estimbounds 
  estimates_girls_nm[[subject]]<-GetBounds(leebounds_wout_monotonicity(leedata_girls,p.0.star[mydata$SEX_NAME==0]))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,mydata=leedata_girls,N_rep=N_rep,p.0.star[mydata$SEX_NAME==0])
  ## confidence region for identified set
  CR_girls_nm[[subject]]<-compute_confidence_region(bounds_bb,estimates_girls_nm[[subject]], ci_alpha=ci_alpha )
}
#p0_girls<-mean(leedata_girls$selection[leedata_girls$treat==0])/mean(leedata_girls$selection[leedata_girls$treat==1])
report_results(estimates_boys_nm,CR_boys_nm)
report_results(estimates_girls_nm,CR_girls_nm)


########################## COLUMN 4: ORTHOGONAL ESTIMATES WITH SIMPLE LOGISTIC REGRESSION#####################
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  
  orthoestimates[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,s.hat=s.hat,sample_size=sample_size)
  CR_ortho[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates[[subject]],ci_alpha=ci_alpha)
  
  
}



for (subject in c("MATH","READING","WRITING")) {
  
  
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  
  orthoestimates[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  
  
  
  leebounds_ortho_result_boys<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==1,],
                                                                 sample_size=sum(leedata_cov$SEX_NAME==1),s.hat=s.hat[leedata_cov$SEX_NAME==1,],
                                                                 quantile_grid_size = quantile_grid_size,
                                                                 y.hat=y.hat[leedata_cov$SEX_NAME==1,])
  
  orthoestimates_boys[[subject]]<-sort(GetBounds(leebounds_ortho_result_boys))
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==1,],B=N_rep,
                                        function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==1,],s.hat=s.hat[leedata_cov$SEX_NAME==1,],sample_size=sum(leedata_cov$SEX_NAME==1))
  CR_ortho_boys[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_boys[[subject]],ci_alpha=ci_alpha)
  
  
  
  
  leebounds_ortho_result_girls<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==0,],
                                                                  sample_size=sum(leedata_cov$SEX_NAME==0),s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                                                                  quantile_grid_size = quantile_grid_size,
                                                                  variables_for_outcome=c("AGE2","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"),
                                                                  y.hat=y.hat[leedata_cov$SEX_NAME==0,])
  
  orthoestimates_girls[[subject]]<-sort(GetBounds(leebounds_ortho_result_girls))
  y.hat<-leebounds_ortho_result$y.hat
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==0,],B=N_rep,
                                        function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==0,],s.hat=s.hat[leedata_cov$SEX_NAME==0,],sample_size=sum(leedata_cov$SEX_NAME==0))
  
  CR_ortho_girls[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_girls[[subject]],ci_alpha=ci_alpha)
  
}
report_results(orthoestimates_boys, CR_ortho_boys)
report_results(orthoestimates_girls, CR_ortho_girls)



############################ COLUMN 7: RESULTS FOR COVARIATE GROUPS WITH POSITIVE LOWER BOUND ########################################



leedata<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK)
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  selected_vars<-c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE","treat","selection","outcome")
  estimates_plb[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  
  res<-summary_subjects_positive_lower_bound(leedata_cov=leedata_cov[,selected_vars],y.hat=y.hat,s.hat=s.hat,inds_helps=(p.0.star<=1)) 
  frac_positive[[subject]]<-res$fraction
  estimates_plb[[subject]]<-GetBounds(res)
  
  res_bb<-weighted_bb(mydata=leedata_cov[,selected_vars],function_name=summary_subjects_positive_lower_bound,
                      y.hat=y.hat,s.hat=s.hat,inds_helps=(p.0.star<=1),B=N_rep)
  
  res_bb<-res_bb[, !is.na(apply(res_bb,2,sum))]
  CR_plb[[subject]]<-compute_confidence_region(ATE_boot=t(res_bb),ATE_est=  estimates_plb[[subject]],ci_alpha=ci_alpha)
  
  
  
  res_boys<-summary_subjects_positive_lower_bound(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==1,selected_vars],
                                                  y.hat=y.hat[leedata_cov$SEX_NAME==1,],
                                                  s.hat=s.hat[leedata_cov$SEX_NAME==1,],
                                                  inds_helps=(p.0.star[leedata_cov$SEX_NAME==1]<=1)) 
  
  estimates_plb_boys[[subject]]<-sort(GetBounds(res_boys))
  
  res_bb_boys<-weighted_bb(mydata=leedata_cov[leedata_cov$SEX_NAME==1,selected_vars],
                           function_name=summary_subjects_positive_lower_bound,
                           y.hat=y.hat[leedata_cov$SEX_NAME==1,],
                           s.hat=s.hat[leedata_cov$SEX_NAME==1,],inds_helps=(p.0.star[leedata_cov$SEX_NAME==1]<=1),B=N_rep)
  res_bb_boys<-res_bb_boys[,!is.na( apply(res_bb_boys,2,sum))]
  CR_plb_boys[[subject]]<-compute_confidence_region(ATE_boot=t(res_bb_boys),ATE_est=  estimates_plb_boys[[subject]],ci_alpha=ci_alpha)
  
  
  
  res_girls<-summary_subjects_positive_lower_bound(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==0,selected_vars],y.hat=y.hat[leedata_cov$SEX_NAME==0,],
                                                   s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                                                   inds_helps=(p.0.star[leedata_cov$SEX_NAME==0]<=1)) 
  
  
  estimates_plb_girls[[subject]]<-sort(GetBounds(res_girls))
  
  res_bb_girls<-weighted_bb(mydata=leedata_cov[leedata_cov$SEX_NAME==0,selected_vars],
                            function_name=summary_subjects_positive_lower_bound,
                            y.hat=y.hat[leedata_cov$SEX_NAME==0,],
                            s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                            inds_helps=(p.0.star[leedata_cov$SEX_NAME==0]<=1),B=N_rep)
  
  res_bb_girls<-res_bb_girls[,!is.na( apply(res_bb_girls,2,sum))]
  CR_plb_girls[[subject]]<-compute_confidence_region(ATE_boot=t(res_bb_girls),ATE_est=  estimates_plb_girls[[subject]],ci_alpha=ci_alpha)
  
}


############################ COLUMN 5: POST-LASSO RESULTS ########################################
form_nonmonotone<-as.formula("selection ~ (treat + AGE2 + SEX_NAME +MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE) * (AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                             DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                             STRATA3 + STRATA4 )")
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])
glm.fit<-estimate_selection(form=form_nonmonotone,selection_function = rlassologit,
                            leedata=leedata_cov,variables_for_selection = setdiff(colnames(leedata_cov),"outcome"))
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
mean(p.0.star<=1)
hist(p.0.star)
glm.fit$coefficients[glm.fit$coefficients!=0]
summary(glm.fit)







for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  
  orthoestimates_postlasso[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,
                                        s.hat=s.hat,sample_size=sample_size)
  
  CR_ortho_postlasso[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_postlasso[[subject]],ci_alpha=ci_alpha)
  
  
  leebounds_ortho_result_boys<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==1,],
                                                                 sample_size=sum(leedata_cov$SEX_NAME==1),
                                                                 s.hat=s.hat[leedata_cov$SEX_NAME==1,],
                                                                 quantile_grid_size = quantile_grid_size,
                                                                 variables_for_outcome=c("AGE2","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  orthoestimates_postlasso_boys[[subject]]<-sort(GetBounds( leebounds_ortho_result_boys))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==1,],
                                        B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==1,],
                                        s.hat=s.hat[leedata_cov$SEX_NAME==1,],sample_size=sum(leedata_cov$SEX_NAME==1))
  
  CR_ortho_postlasso_boys[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_postlasso_boys[[subject]],ci_alpha=ci_alpha)
  
  leebounds_ortho_result_girls<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==0,],
                                                                  sample_size=sum(leedata_cov$SEX_NAME==0),
                                                                  s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                                                                  quantile_grid_size = quantile_grid_size,
                                                                  variables_for_outcome=c("AGE2","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  orthoestimates_postlasso_girls[[subject]]<-sort(GetBounds( leebounds_ortho_result_girls))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==0,],
                                        B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==0,],
                                        s.hat=s.hat[leedata_cov$SEX_NAME==0,],sample_size=sum(leedata_cov$SEX_NAME==0))
  
  CR_ortho_postlasso_girls[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_postlasso_girls[[subject]],ci_alpha=ci_alpha)
  
}



##################### RANDOM FOREST: COLUMN 6 ##############################################

sample_size<-dim(mydata)[1]
inds_train<-sample(1:sample_size,floor(dim(mydata)[1]/2))
inds_test<-setdiff(1:sample_size,inds_train)
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])
leedata_cov$selection<-as.logical(leedata_cov$selection)
form_rf<-as.formula(paste0("selection~(", paste0(exogenous_covariates,collapse="+"),")+treat"))
rf.fit.train<-ranger(form_rf,leedata_cov[inds_train,],importance = "permutation",num.trees=30,probability=TRUE)
s.hat<-data.frame(s.0.hat=rep(NA,sample_size),s.1.hat=rep(NA,sample_size))
s.hat[inds_test,]<-predict_selection(rf.fit.train,leedata_cov[inds_test,c("treat","selection",exogenous_covariates)],rf=TRUE)


rf.fit.test<-ranger(form_rf,leedata_cov[inds_test,],importance = "permutation",num.trees=30,probability=TRUE)
s.hat[inds_train,]<-predict_selection(rf.fit.test,leedata_cov[inds_train,c("treat","selection",exogenous_covariates)],rf=TRUE)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat


for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME","MOM_SCH","DAD_SCH","DAD_AGE","MOM_AGE"))
  
  orthoestimates_rf[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,s.hat=s.hat,sample_size=sample_size)
  CR_ortho_rf[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf[[subject]],ci_alpha=ci_alpha)
  
  
  
  leebounds_ortho_result_boys<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==1,],
                                                                 sample_size=sum(leedata_cov$SEX_NAME==1),s.hat=s.hat[leedata_cov$SEX_NAME==1,],
                                                                 y.hat=y.hat[leedata_cov$SEX_NAME==1,])
  
  orthoestimates_rf_boys[[subject]]<-sort(GetBounds( leebounds_ortho_result_boys))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==1,],B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==1,],s.hat=s.hat[leedata_cov$SEX_NAME==1,],sample_size=sum(leedata_cov$SEX_NAME==1))
  CR_ortho_rf_boys[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf_boys[[subject]],ci_alpha=ci_alpha)
  
  
  leebounds_ortho_result_girls<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==0,],
                                                                  sample_size=sum(leedata_cov$SEX_NAME==0),s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                                                                  y.hat=y.hat[leedata_cov$SEX_NAME==0,])
  
  orthoestimates_rf_girls[[subject]]<-sort(GetBounds( leebounds_ortho_result_girls))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==0,],B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==0,],s.hat=s.hat[leedata_cov$SEX_NAME==0,],sample_size=sum(leedata_cov$SEX_NAME==0))
  CR_ortho_rf_girls[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf_girls[[subject]],ci_alpha=ci_alpha)
}



#Save results
save.image(paste0(my_path,"/ABBKK_2002/data/ABBKK.RData"))
