## clear workspace
rm(list=ls())
## load relevant libraries
library(expm)
library(hdm)
## specify path to code 
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/aux.R"))
## load data
mydata<-read.csv(paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
## clean data: take those wit age <=25 (age2<=23) and non-missing record of sex name
mydata<-mydata[mydata$AGE2<=23 & mydata$SEX_IS_NA == 0, ]
dim(mydata)

#Parameters
sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002.txt"),append=TRUE)

## confidence level
ci_alpha=0.05
## number of bootstrap repetition
N_rep=1000
## grid size for quantile regression
quantile_grid_size=0.01
# Save results here

frac_positive<-list()
estimates_plb<-list()
CR_plb<-list()

orthoestimates_postlasso<-list()
CR_ortho_postlasso<-list()


print ("Parameters:")
print(paste0("Confidence level ",ci_alpha))
print(paste0("Number of bootstrap repetitions ",N_rep))
print(paste0("Quantile grid size for quantile regression ",quantile_grid_size))
## exogenous covariates from ABBKK data
exogenous_covariates<-c("AGE2","SEX_NAME","MOM_SCH" ,"MOM_AGE","DAD_SCH" ,"DAD_AGE",
                        "DAREA4",   "DAREA5",   "DAREA6",   "DAREA7",
                        "DAREA11" ,  "DAREA15",
                        "DAREA16" , "DAREA17"  ,"DAREA18" , "DAREA19",
                        
                        "STRATA1",  "STRATA2",  "STRATA3" , "STRATA4", "MOM_AGE_IS_NA","DAD_AGE_IS_NA","MOM_SCH_IS_NA","DAD_SCH_IS_NA"
)
print("Estimating test participation by post-lasso-logistic")
form_nonmonotone<-as.formula("selection ~ (treat + AGE2 + SEX_NAME +MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE) * (AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                             DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                             STRATA3 + STRATA4 +MOM_AGE_IS_NA+DAD_AGE_IS_NA)")
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])
glm.fit<-estimate_selection(form=form_nonmonotone,selection_function = rlassologit,
                            leedata=leedata_cov,variables_for_selection = setdiff(colnames(leedata_cov),"outcome"))
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
s.hat.lasso<-s.hat

## split interactions into raw covariates
vars_for_selection<-setdiff(unique( unlist( strsplit (c(names(glm.fit$coefficients[glm.fit$coefficients!=0])
),split=":"))), c("(Intercept)", "treat" ) )

print (paste0("Fraction of subjects whose conditional trimming threshold is less than one (post-lasso-logistic) ", round (mean(p.0.star<=1),3)))
print(paste0("Selected covariates for test participation are ",paste0(setdiff(vars_for_selection,c("(Intercept)","treat") ),collapse=" ")))

print("Selecting covariates for test score by post-lasso")
form_outcome<-as.formula("outcome ~ treat + (AGE2 + SEX_NAME +MOM_SCH + MOM_AGE + DAD_SCH + 
                         DAD_AGE) * (AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                         DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                         DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                         STRATA3 + STRATA4 +MOM_AGE_IS_NA+DAD_AGE_IS_NA)")
vars_for_outcome<-c()
leedata<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK)
for (subject in c("MATH","READING","WRITING")) {
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  fit<-rlasso(form_outcome,leedata_cov[leedata_cov$selection==1 & leedata_cov$treat==0,])
  vars_for_outcome<-c(vars_for_outcome,names(fit$coefficients[fit$coefficients!=0]))
}
vars_for_outcome<-setdiff(vars_for_outcome,"(Intercept)")

## split interactions into raw covariates
vars_for_outcome<-setdiff(unique( unlist( strsplit (c(vars_for_outcome
),split=":"))), c("(Intercept)", "treat" ) )

print(paste0("Selected covariates for test score are ",paste0(vars_for_outcome,collapse=" ")))
print(paste0("Manually added covariates for test score are ",c("AGE2","MOM_AGE","DAD_AGE_IS_NA"),collapse=" "))



vars_for_outcome<-c("AGE2","SEX_NAME","MOM_AGE","DAD_AGE_IS_NA",
                    "MOM_AGE_IS_NA", "DAD_SCH","MOM_SCH","DAD_AGE","DAREA4","DAREA6","DAREA15")
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,s.hat=s.hat.lasso,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=vars_for_outcome)
  
  orthoestimates_postlasso[[subject]]<-GetBounds(leebounds_ortho_result)
  y.hat<-leebounds_ortho_result$y.hat

  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,
                                        s.hat=s.hat.lasso)
  
  CR_ortho_postlasso[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_postlasso[[subject]],ci_alpha=ci_alpha)
  
} 


############################ COLUMN 6: RESULTS FOR COVARIATE GROUPS WITH POSITIVE LOWER BOUND ########################################


for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,
                                                            s.hat=s.hat.lasso,
                                                            quantile_grid_size = quantile_grid_size,variables_for_outcome = exogenous_covariates)
  y.hat<-leebounds_ortho_result$y.hat
  
  res<-summary_subjects_positive_lower_bound(leedata_cov=leedata_cov,
                                             y.hat=y.hat,
                                             s.hat=s.hat.lasso,
                                             inds_helps=(p.0.star<=1)) 
  frac_positive[[subject]]<-res$fraction
  estimates_plb[[subject]]<-GetBounds(res)
  
  res_bb<-weighted_bb(mydata=leedata_cov,function_name=summary_subjects_positive_lower_bound,
                      y.hat=y.hat,s.hat=s.hat.lasso,inds_helps=(p.0.star<=1),B=N_rep)
  
  res_bb<-res_bb[, !is.na(apply(res_bb,2,sum)) ]
  CR_plb[[subject]]<-compute_confidence_region(ATE_boot=t(res_bb),ATE_est=  estimates_plb[[subject]],ci_alpha=ci_alpha)
  
}

## ABBKK 2002, Table 1 column 5
print("Column 5")
report_results(orthoestimates_postlasso,CR_ortho_postlasso)
## ABBKK 2002, Table 1 column 6
print("Column 6")
report_results(estimates_plb,CR_plb)
