## clear workspace
rm(list=ls())
## load relevant libraries
library(expm)
## specify path to code 
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/auxiliary.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/orthogonal_correction.R"))
source(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/aux.R"))
## load data
mydata<-read.csv(paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
## open log file
sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002_baseline.txt"),append=FALSE)


## clean data: take those wit age <=25 (age2<=23) and non-missing record of sex name
mydata<-mydata[mydata$AGE2<=23 & mydata$SEX_IS_NA == 0 , ]
dim(mydata)
#Parameters


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
                        
                        "STRATA1",  "STRATA2",  "STRATA3" , "STRATA4", "DAD_AGE_IS_NA","MOM_SCH_IS_NA","MOM_AGE_IS_NA","DAD_SCH_IS_NA"
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

orthoestimates2<-list()
CR_ortho2<-list()

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
form_nonmonotone<-as.formula(paste0("selection ~ treat *(",paste0(exogenous_covariates,collapse="+"), ")"))

leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])


glm.fit<-glm(form=form_nonmonotone,
             data=leedata_cov[, setdiff(colnames(leedata_cov),"outcome")],
             family="binomial")

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
  bounds_bb<-bounds_bb[!is.na(apply(bounds_bb,1,sum)),]
  ## confidence region for identified set
  CR_nonmonotone[[subject]]<-compute_confidence_region(bounds_bb,estimates_nonmonotone[[subject]], ci_alpha=ci_alpha )
}

########################## COLUMN 4: ORTHOGONAL ESTIMATES WITH SIMPLE LOGISTIC REGRESSION (POWER 1)#####################

print ("Estimating Column 4: Lee (2009) bounds without monotonicity and 24 covariates ... ")
for (subject in c("MATH","READING","WRITING")) {
  
  
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  # estimate quantile regression and assemble conditional quantile
  first_stage<-first_stage_wrapper(leedata_cov=leedata_cov,s.hat=s.hat.logistic,
                                   quantile_grid_size = quantile_grid_size,
                                   variables_for_outcome=exogenous_covariates)
  
  
  leebounds_ortho_result<-second_stage_wrapper(leedata=leedata_cov,
                                               inds_helps=first_stage$inds_helps,
                                               y.hat=first_stage$y.hat,
                                               s.hat=s.hat,ortho=FALSE)
  
  orthoestimates[[subject]]<-GetBounds(leebounds_ortho_result)
  

  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=second_stage_wrapper,
                                        y.hat=first_stage$y.hat,
                                        inds_helps=first_stage$inds_helps,
                                        s.hat=s.hat.logistic,ortho=TRUE,c_quant=1)
  
  CR_ortho[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates[[subject]],ci_alpha=ci_alpha)
  
}

########################## COLUMN 5: ORTHOGONAL ESTIMATES WITH SIMPLE LOGISTIC REGRESSION (POWER 2)#####################

form_nonmonotone<-as.formula(paste0("selection ~ (treat+AGE2  +MOM_SCH + MOM_AGE + DAD_SCH + 
                                    DAD_AGE) *(",paste0(exogenous_covariates,collapse="+"), ")" ))

leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])


glm.fit<-glm(form=form_nonmonotone,
             data=leedata_cov[, setdiff(colnames(leedata_cov),"outcome")],
             family="binomial")
s.hat<-predict_selection(glm.fit,leedata_cov[,c("treat","selection",exogenous_covariates)])
s.hat<-data.frame(s.0.hat=s.hat$s.0.hat, s.1.hat=s.hat$s.1.hat)
s.hat.logistic<-s.hat

p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
print (paste0("Fraction of subjects whose conditional trimming threshold is less than one (logistic regression)", round (mean(p.0.star<=1),3)))


print ("Estimating Column 5: Lee (2009) bounds without monotonicity and 24*5 covariates ... ")
for (subject in c("MATH","READING","WRITING")) {
  
  
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  # estimate quantile regression and assemble conditional quantile
  first_stage<-first_stage_wrapper(leedata_cov=leedata_cov,s.hat=s.hat.logistic,
                                   quantile_grid_size = quantile_grid_size,
                                   variables_for_outcome=exogenous_covariates)
  
  
  leebounds_ortho_result<-second_stage_wrapper(leedata=leedata_cov,
                                               inds_helps=first_stage$inds_helps,
                                               y.hat=first_stage$y.hat,
                                               s.hat=s.hat,ortho=FALSE)
  
  orthoestimates2[[subject]]<-GetBounds(leebounds_ortho_result)

  
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,
                                        B=N_rep,function_name=second_stage_wrapper,
                                        y.hat=first_stage$y.hat,
                                        inds_helps=first_stage$inds_helps,
                                        s.hat=s.hat.logistic,ortho=TRUE,c_quant=1)
  
  CR_ortho2[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates2[[subject]],ci_alpha=ci_alpha)
  
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
## ABBKK 2002, Table 1 column 5
print("Column 5")
report_results(orthoestimates2,CR_ortho2)

## stop sinking
sink(file=NULL)
closeAllConnections()

### width computation
compute_width<-function(x) {
  z=x[2]-x[1]
  return(round(z,3))
}
width_estimates_2=unlist(lapply(estimates,compute_width))
width_CR_2=unlist(lapply(CR,compute_width))

width_estimates_3=unlist(lapply(estimates_nonmonotone,compute_width))
width_CR_3=unlist(lapply(CR_nonmonotone,compute_width))

width_estimates_4=unlist(lapply(orthoestimates,compute_width))
width_CR_4=unlist(lapply(CR_ortho,compute_width))

width_estimates_5=unlist(lapply(orthoestimates2,compute_width))
width_CR_5=unlist(lapply(CR_ortho2,compute_width))

sink(paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table_ABBKK_2002.txt"),append=TRUE)
print("Width computation: Column 2")
width_estimates_2
width_CR_2
print("Width computation: Column 3")
width_estimates_3
width_CR_3
print("Width computation: Column 4")
width_estimates_4
width_CR_4
print("Width computation: Column 5")
width_estimates_5
width_CR_5

save.image(paste0(my_path,"/ABBKK_2002/draft_older/ABBKK_2002_baseline.RData"))