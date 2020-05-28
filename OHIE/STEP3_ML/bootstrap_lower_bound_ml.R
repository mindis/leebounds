## clear workspace

## this code replicates Col 1-2, 4-5, 7-8 of Table Utilization and Health, Semenova (2020)
rm(list=ls())


library(feather)
library(tidyverse)
## for weighted standard errors
library(radiant.data)
## for cluster-robust ITT
library(miceadds)
## for LATE
library(AER)
library(ivpack)
library(hdm)

## for parallel computing
library(doParallel)
library(foreach)
cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))

my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path, "/R/auxiliary.R"))
source(paste0(my_path, "/R/leebounds.R"))

source(paste0(my_path, "/OHIE/STEP3_ML/ml_utils.R"))
########### READ IN DATA ################

prepared_data<-read_feather(paste0(my_path,"/OHIE/OHIE_Data/prepared_data.feather"))
prepared_data_sample_12m<-filter(prepared_data,sample_12m=="In 12m mail survey sample")
prepared_data_responded_control<-filter(prepared_data_sample_12m,treatment == 0)

stratification_controls<-grep("ddd",colnames(prepared_data),value=TRUE)
weight_name<-"weight_12m"
############# OUTCOMES #########


discrete_outcome_names<-c("rx_any_12m","doc_any_12m","er_any_12m","hosp_any_12m",
                          "cost_any_oop_12m","cost_any_owe_12m","cost_borrow_12m","cost_refused_12m", 
                          "health_genflip_bin_12m","health_notpoor_12m","health_chgflip_bin_12m","nodep_screen_12m")


trimming_value_list<-c(0,0,0,0,1,1,1,1,0,0,0,0)

utilization_outcomes<-c("rx_num_mod_12m","doc_num_mod_12m","er_num_mod_12m","hosp_num_mod_12m")
cost_outcomes<-c(   "cost_doc_oop_12m",        "cost_er_oop_12m"  ,       "cost_rx_oop_12m",         "cost_oth_oop_12m" ,       "cost_tot_oop_12m" ,      
                    "cost_tot_owe_12m" ,          "cost_tot_oop_mod_12m" ,   "cost_tot_owe_mod_12m")

numeric_health_outcomes<-c("notbaddays_phys_12m","notbaddays_tot_12m","notbaddays_ment_12m")
all_names<-c(utilization_outcomes,cost_outcomes,numeric_health_outcomes,discrete_outcome_names)

#,"group_id"
#### all covariates = compulsory (stratification) and additional (pre-randomization) covariates
list_of_all_controls<-setdiff(colnames(prepared_data), c("person_id","household_id","treatment", "ohp_all_ever_survey",
                                                         all_names, "weight_12m","sample_12m","sample_12m_resp","group_id"
))

### compulsory covariates
stratification_controls<-grep("ddd",colnames(prepared_data_sample_12m),value=TRUE)
### covariates for point mass estimation
selected_covariates_all_regressions<-c("female_list","english_list","zip_msa","num_visit_pre_cens_ed",
                                       "num_out_pre_cens_ed","snap_ever_prenotify07","snap_tot_hh_prenotify07","tanf_ever_prenotify07","tanf_tot_hh_prenotify07",
                                       stratification_controls[c(1:7,9,11,13:16)])


weight_name<-"weight_12m"
B<-400
############# Estimate bounds ############
trimming_value_list<-c(0,0,0,0,1,1,1,1,0,0,0,0)
quantile_grid_size=0.01

cluster_robust_se<-list()
B=400

formula_threshold_list<-list()
formula_point_mass_list<-list()

formula_threshold_list[["discrete"]]<-as.formula(paste0("selection~treat*(female_list
                                                        +zip_msa+snap_tot_hh_prenotify07+
                                                        ed_charg_tot_pre_ed+snap_ever_prenotify07+
                                                        charg_tot_pre_ed+num_epct_pre_ed+
                                                        num_out_pre_cens_ed+num_visit_pre_cens_ed+
                                                        any_out_pre_ed+num_edcnnp_pre_ed+
                                                        any_hiun_pre_ed+num_hiun_pre_cens_ed+
                                                        num_ne_pre_ed+num_on_pre_cens_ed+
                                                        tanf_tot_hh_prenotify07+any_on_pre_ed+
                                                        any_visit_pre_ed)*(",paste0(list_of_all_controls,collapse="+"),")" ))

formula_point_mass_list[["discrete"]]<-as.formula(paste0(paste0("outcome_point_mass~",paste0(c("female_list","english_list","zip_msa","num_visit_pre_cens_ed",
                                                                            "num_out_pre_cens_ed","snap_ever_prenotify07","snap_tot_hh_prenotify07","tanf_ever_prenotify07","tanf_tot_hh_prenotify07",
                                                                            stratification_controls),collapse="+") )))




formula_threshold_list[["utilization"]]<-as.formula(paste0("selection~treat+ ",paste0(selected_covariates_all_regressions,collapse="+") ))
formula_point_mass_list[["utilization"]]<-as.formula(paste0("outcome_point_mass~ ",paste0(selected_covariates_all_regressions,collapse="+") ))

formula_threshold_list[["numeric_health"]]<-as.formula(paste0("selection~treat+ ",paste0(list_of_all_controls,collapse="+") ))
formula_point_mass_list[["numeric_health"]]<-as.formula(paste0("outcome_point_mass~ ",paste0(stratification_controls,collapse="+") ))

prepared_data_sample_12m<-as.data.frame(prepared_data_sample_12m)
prepared_data_sample_12m<-sapply(prepared_data_sample_12m,as.numeric)
prepared_data_sample_12m<-as.data.frame(prepared_data_sample_12m)
cluster_robust_se<-list()
B=400
#myres<-matrix(NA,3,B)
for (outcome_name in utilization_outcomes) {
  
  print(outcome_name)
 
  
  myres<- foreach(b=1:B, .combine = 'cbind',.packages = c("tidyverse","AER","miceadds","hdm"))  %dopar% {
    ## some (a small number) bootstrap draws may fail the basic trimming method
    #for (b in 1:B) {
     ans<-bootstrap_mean_itt_late_ml(b,
                                    first_stage_wrapper_function = first_stage_wrapper_continuous,
                                    outcome_name=outcome_name,
                                    list_of_all_controls=list_of_all_controls,
                                    selection_function_for_threshold=glm,
                                    selection_function_for_threshold_name="glm",
                                    stratification_controls= stratification_controls[c(1:7,9,11)],
                                    selection_function_point_mass=rlassologit,
                                    selection_function_point_mass_name="rlassologit",
                                    formula_threshold=formula_threshold_list[["utilization"]],
                                    formula_point_mass=formula_point_mass_list[["utilization"]],
                                    variables_for_outcome=setdiff(selected_covariates_all_regressions,stratification_controls[c(8,10,12:16)]),
                                    quantile_grid_size=quantile_grid_size,
                                    sort_quantiles=TRUE,
                                    selected_inds=NULL,
                                    max_p0_hat=1.001001,
                                    point_mass_cutoff=0.2,
                                    point_mass=TRUE,
                                    weight_name="weight_12m",
                                    trimming_value=0)
    #
  }
  cluster_robust_se[[outcome_name]]<-apply(myres,1,sd)
 }

for (outcome_name in numeric_health_outcomes) {
  
  print(outcome_name)
  myres<- foreach(b=1:B, .combine = 'cbind',.packages = c("tidyverse","AER","miceadds","hdm"))  %dopar% {
    ## some (a small number) bootstrap draws may fail the basic trimming method
    ans<-tryCatch(bootstrap_mean_itt_late_ml(b,
                                    first_stage_wrapper_function = first_stage_wrapper_continuous,
                                    outcome_name=outcome_name,
                                    list_of_all_controls=list_of_all_controls,
                                    selection_function_for_threshold=glm,
                                    selection_function_for_threshold_name="glm",
                                    stratification_controls= stratification_controls[c(1:7,9,11)],
                                    selection_function_point_mass=rlassologit,
                                    selection_function_point_mass_name="rlassologit",
                                    formula_threshold=formula_threshold_list[["numeric_health"]],
                                    formula_point_mass=formula_point_mass_list[["numeric_health"]],
                                    variables_for_outcome=stratification_controls[c(1:7,9,11)],
                                    quantile_grid_size=0.01,
                                    sort_quantiles=TRUE,
                                    selected_inds=NULL,
                                    max_p0_hat=1.001001,
                                    point_mass_cutoff=0.20,
                                    point_mass=TRUE, trimming_value=0,weight_name="weight_12m"),error=function(e) NULL)
    #tryCatch(
  }
  cluster_robust_se[[outcome_name]]<-apply(myres,1,sd)
  
}


for (outcome_name in discrete_outcome_names) {
    
  print(outcome_name)
  myres<- foreach(b=1:B, .combine = 'cbind',.packages = c("tidyverse","AER","miceadds","hdm"))  %dopar% {
    ## some (a small number) bootstrap draws may fail the basic trimming method
    ans<-bootstrap_mean_itt_late_ml(b,
                           
                                             first_stage_wrapper_function = first_stage_wrapper_discrete,
                                             outcome_name=outcome_name,
                                             list_of_all_controls=list_of_all_controls,
                                             trimming_value=trimming_value_list[match(outcome_name,discrete_outcome_names)],
                                             selection_function_for_threshold=rlassologit,
                                             selection_function_for_threshold_name="rlassologit",
                                             stratification_controls= stratification_controls[c(1:7,9,11)],
                                             selection_function_point_mass=rlassologit,
                                             selection_function_point_mass_name="rlassologit",
                                             formula_threshold=formula_threshold_list[["discrete"]],
                                             formula_point_mass=formula_point_mass_list[["discrete"]],
                                             weight_name="weight_12m",point_mass_cutoff=0.20,
                                             max_p0_hat=1.000001)
    #tryCatch(,error=function(e) NULL)
  }

  cluster_robust_se[[outcome_name]]<-apply(myres,1,sd)
 }



ml_trim_sd<-matrix(unlist(cluster_robust_se),nrow=3,byrow=FALSE)
rownames(ml_trim_sd)<-c("control_mean","itt","late")
#colnames(ml_trim_sd)<-c(utilization_outcomes,numeric_health_outcomes,discrete_outcome_names[1:4],discrete_outcome_names[9:12])
colnames(ml_trim_sd)<-utilization_outcomes

ml_trim_sd<-t(ml_trim_sd)
ml_trim_sd<-apply(ml_trim_sd,2,round,4)



write.csv(ml_trim_sd,paste0(my_path,"/OHIE/STEP3_ML/csv/sd_ml_trim_utilization.csv"))
#save.image(paste0(my_path,"/OHIE/STEP3_ML/draft/ml_trim_boot_utilization.RData"))
