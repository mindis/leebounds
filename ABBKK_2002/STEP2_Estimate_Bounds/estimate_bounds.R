## clear workspace
rm(list=ls())
## load relevant libraries
library(expm)
## specify path to code 
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/aux.R"))
## load data
mydata<-read.csv(paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
dim(mydata)
## clean data: take those wit age <=25 (age2<=23) and non-missing record of sex name
mydata<-mydata[mydata$AGE2<=23 & mydata$SEX_IS_NA == 0, ]

#Parameters
sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002.txt"),append=FALSE)

## confidence level
ci_alpha=0.05
## number of bootstrap repetition
N_rep=1000
## grid size for quantile regression
quantile_grid_size=0.01

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

## Save Results Column 1
## estimates
estimates<-list()
## confidence region
CR<-list()

estimates_nonmonotone<-list()
CR_nonmonotone<-list()

orthoestimates<-list()
CR_ortho<-list()

leedata<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK)


############### COLUMN 2: BASIC LEE (2008) BOUNDS UNDER MONOTONICITY ####################
print ("Estimating Column 2: Lee (2009) bounds under monotonicity with no covariates ... ")
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
print(paste0(" Unconditional trimming threshold ",round(p0,3)))
################ COLUMN 3: LEE (2008) bounds without mononotonicity and covariates ####################

## 
form_nonmonotone<-as.formula("selection ~ treat * (AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                             DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                             DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                             STRATA3 + STRATA4 +MOM_AGE_IS_NA+DAD_AGE_IS_NA+MOM_SCH_IS_NA+DAD_SCH_IS_NA )")
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])


glm.fit<-estimate_selection(form=form_nonmonotone,selection_function = glm,
                            leedata=leedata_cov,variables_for_selection = setdiff(colnames(leedata_cov),"outcome"))
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
s.hat.logistic<-s.hat

p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
print (paste0("Fraction of subjects whose conditional trimming threshold is less than one (logistic regression)", round (mean(p.0.star<=1),3)))


print ("Estimating Column 3: Lee (2009) bounds without monotonicity and covariates ... ")
for (subject in c("MATH","READING","WRITING")) {
  leedata$outcome<-mydata[,subject]
  
  ## compute estimbounds 
  estimates_nonmonotone[[subject]]<-GetBounds(leebounds_wout_monotonicity(leedata,p.0.star=p.0.star))
  ##  bootstrap draws of bounds
  bounds_bb<-main_bb(function_name=leebounds_wout_monotonicity,mydata=leedata,N_rep=N_rep,p.0.star=p.0.star)
  ## confidence region for identified set
  CR_nonmonotone[[subject]]<-compute_confidence_region(bounds_bb,estimates_nonmonotone[[subject]], ci_alpha=ci_alpha )
}

########################## COLUMN 4: ORTHOGONAL ESTIMATES WITH SIMPLE LOGISTIC REGRESSION#####################

print ("Estimating Column 4: Lee (2009) bounds without monotonicity and 25 covariates ... ")
for (subject in c("MATH","READING","WRITING")) {
  
  
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,s.hat=s.hat.logistic,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=exogenous_covariates)
  
  orthoestimates[[subject]]<-GetBounds(leebounds_ortho_result)
  y.hat<-leebounds_ortho_result$y.hat
  
  
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,
                                        s.hat=s.hat.logistic)
  
  CR_ortho[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates[[subject]],ci_alpha=ci_alpha)
  
}
##### ###   REPORT RESULTS #############


## ABBKK 2002, Table 1 column 2
print("Column 2")
report_results(estimates,CR)
## ABBKK 2002, Table 1 column 3
print("Column 3")
report_results(estimates_nonmonotone,CR_nonmonotone)
## ABBKK 2002, Table 1 column 4
print("Column 4")
report_results(orthoestimates,CR_ortho)
