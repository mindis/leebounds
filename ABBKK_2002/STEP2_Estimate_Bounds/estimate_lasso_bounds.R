## clear workspace
rm(list=ls())
## load relevant libraries
library(expm)
library(hdm)
## specify path to code 
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/orthogonal_correction.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/aux.R"))
## load data
mydata<-read.csv(paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
## clean data: take those wit age <=25 (age2<=23) and non-missing record of sex name
mydata<-mydata[mydata$AGE2<=23 & mydata$SEX_IS_NA == 0, ]
dim(mydata)

#Parameters
#sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002_baseline.txt"),append=TRUE)

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
                        
                        "STRATA1",  "STRATA2",  "STRATA3" , "STRATA4","MOM_AGE_IS_NA","DAD_AGE_IS_NA","MOM_SCH_IS_NA","DAD_SCH_IS_NA"
)
print("Estimating test participation by post-lasso-logistic")
## preliminary  selection of geographical indicators
form_nonmonotone<-as.formula(paste0("selection ~ (treat+AGE2  +MOM_SCH + MOM_AGE + DAD_SCH + 
                                    DAD_AGE) * (AGE2  +MOM_SCH + MOM_AGE + DAD_SCH + 
                                    DAD_AGE)*(",paste0(exogenous_covariates,collapse="+"),")"))

leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])
glm.fit<-rlassologit(form=form_nonmonotone,
                     data=leedata_cov[,setdiff(colnames(leedata_cov),"outcome")])
vars_for_selection<-setdiff(unique( unlist( strsplit (c(names(glm.fit$coefficients[glm.fit$coefficients!=0])
),split=":"))), c("(Intercept)", "treat" ) )

print(paste0("Selected geographic indiactors are",paste0(grep("STRATA|DAREA",vars_for_selection,value=TRUE),collapse=" ")))


## two-level specification based on geo indicators selected from above

form_nonmonotone<-as.formula("selection ~ (treat+AGE2 + SEX_NAME +MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE) *(MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE+DAREA11 + DAREA17 + DAREA18+DAREA19 + STRATA2 + 
                             +MOM_AGE_IS_NA+DAD_AGE_IS_NA)")
glm.fit<-estimate_selection(form=form_nonmonotone,leedata=leedata_cov,
                            selection_function=rlassologit,
                            selection_function_name="rlassologit",
                            variables_for_selection =setdiff(colnames(leedata_cov),"outcome") )
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
s.hat.lasso<-s.hat


vars_for_selection<-setdiff(unique( unlist( strsplit (c(names(glm.fit$coefficients[glm.fit$coefficients!=0])
),split=":"))), c("(Intercept)", "treat" ) )
print(paste0("Selected covariates for test participation are ",paste0(setdiff(vars_for_selection,c("(Intercept)","treat") ),collapse=" ")))

print("Selecting covariates for test score by post-lasso")
form_outcome<-as.formula(paste0("outcome ~  (AGE2+MOM_SCH+DAD_SCH+MOM_AGE+DAD_AGE)*(",paste0(exogenous_covariates,collapse="+"), ")" ))
vars_for_outcome<-c()

leedata<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK)

for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  fit<-rlasso(form_outcome,leedata_cov[leedata_cov$selection==1 & leedata_cov$treat==1 ,])
  vars_for_outcome<-c(vars_for_outcome,names(fit$coefficients[fit$coefficients!=0]))
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  fit<-rlasso(form_outcome,leedata_cov[leedata_cov$selection==1 & leedata_cov$treat==0 ,])
  vars_for_outcome<-c(vars_for_outcome,names(fit$coefficients[fit$coefficients!=0]))
}
vars_for_outcome<-setdiff(vars_for_outcome,"(Intercept)")

## split interactions into raw covariates
vars_for_outcome<-setdiff(unique( unlist( strsplit (c(vars_for_outcome
),split=":"))), c("(Intercept)", "treat" ) )

print(paste0("Selected covariates for test score are ",paste0(vars_for_outcome,collapse=" ")))


first_stage_list<-list()
second_stage_list<-list()
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  first_stage<-first_stage_wrapper(leedata_cov=leedata_cov,
                                   s.hat=s.hat.lasso,
                                   quantile_grid_size = quantile_grid_size,
                                   variables_for_outcome=vars_for_outcome)
  first_stage_list[[subject]]<-first_stage
  
  leebounds_ortho_result<-second_stage_wrapper(leedata=leedata_cov,
                                               s.hat=s.hat.lasso,
                                               y.hat=first_stage$y.hat,
                                               inds_helps=first_stage$inds_helps,ortho=TRUE,c_quant=1)
  
  second_stage_list[[subject]]<-leebounds_ortho_result
  orthoestimates_postlasso[[subject]]<-GetBounds(leebounds_ortho_result)
  y.hat<-leebounds_ortho_result$y.hat

  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=second_stage_wrapper,
                                        y.hat=first_stage$y.hat,
                                        inds_helps=first_stage$inds_helps,
                                        s.hat=s.hat.lasso,ortho=TRUE,c_quant=1)
  
  CR_ortho_postlasso[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_postlasso[[subject]],ci_alpha=ci_alpha)
  ### subjects with positive lower bound
  
}
vars_for_plb<-c("MOM_AGE_IS_NA","MOM_AGE","MOM_SCH","MOM_SCH_IS_NA")
CR_plb<-list()
estimates_bb_list<-list()
N_rep=200
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  res<-summary_subjects_positive_lower_bound(leedata_cov[,c("treat","outcome","selection",vars_for_plb)],
                                             form_selection = paste0("selection~(treat)*(",paste0(vars_for_plb,collapse="+"),")"),
                                             selection_function_name="glm",
                                             quantile_grid_size = quantile_grid_size,
                                             variables_for_outcome=vars_for_plb,
                                             form_outcome_plb=as.formula("outcome~."))
  
  estimates_plb[[subject]]<-GetBounds(res)
  frac_positive[[subject]]<-GetFraction(res)
  
  
  
  estimates_bb_list[[subject]]<-weighted_bb(leedata_cov[,c("treat","outcome","selection",vars_for_plb)],
                                            form_selection = paste0("selection~(treat)*(",paste0(vars_for_plb,collapse="+"),")"),
                                            selection_function_name="glm",
                                            quantile_grid_size = quantile_grid_size,
                                            variables_for_outcome=vars_for_plb,
                                            form_outcome_plb=as.formula("outcome~."),function_name = summary_subjects_positive_lower_bound,B=N_rep)
  
  estimates_bb<- estimates_bb_list[[subject]]
  estimates_bb<-estimates_bb[, !is.na(apply(estimates_bb,2,sum))]
  CR_plb[[subject]]<-c(quantile(t(  estimates_bb)[,1],ci_alpha/2,na.rm=TRUE),quantile(t(  estimates_bb)[,2],1-ci_alpha/2,na.rm=TRUE) )
}





## ABBKK 2002, Table 1 column 5
print("Column 6")
report_results(orthoestimates_postlasso,CR_ortho_postlasso)
## ABBKK 2002, Table 1 column 6
print("Column 7")
report_results(estimates_plb,CR_plb)
### width computation
compute_width<-function(x) {
  z=x[2]-x[1]
  return(round(z,3))
}
width_estimates_6=unlist(lapply(orthoestimates_postlasso,compute_width))
width_CR_6=unlist(lapply(CR_ortho_postlasso,compute_width))

width_estimates_7=unlist(lapply(estimates_plb,compute_width))
width_CR_7=unlist(lapply(CR_plb,compute_width))


sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002.txt"),append=TRUE)
print("Width computation: Column 6")
width_estimates_6
width_CR_6
print("Width computation: Column 7")
width_estimates_7
width_CR_7
## stop sinking
sink(file=NULL)
closeAllConnections()
save.image(paste0(my_path,"/ABBKK_2002/draft_older/ABBKK_2002_Col67.RData"))