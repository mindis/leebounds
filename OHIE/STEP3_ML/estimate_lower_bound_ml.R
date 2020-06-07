rm(list=ls())

library(feather)
library(tidyverse)
library(ivpack)
library(miceadds)
library(AER)
## for machine learning part
library(hdm)
########### READ IN DATA ################
my_path<-"/net/holyparkesec/data/tata/leebounds/"

prepared_data<-read_feather(paste0(my_path,"/OHIE/OHIE_Data/prepared_data.feather"))
prepared_data_sample_12m<-filter(prepared_data,sample_12m=="In 12m mail survey sample")

########## READ IN HELPER FUNCTIONS #########
source(paste0(my_path, "/R/auxiliary.R"))
source(paste0(my_path, "/R/leebounds.R"))
source(paste0(my_path, "/OHIE/STEP3_ML/ml_utils.R"))
##### outcome names

discrete_outcome_names<-c("rx_any_12m","doc_any_12m","er_any_12m","hosp_any_12m","cost_any_oop_12m","cost_any_owe_12m","cost_borrow_12m","cost_refused_12m", 
                          "health_genflip_bin_12m","health_notpoor_12m","health_chgflip_bin_12m","nodep_screen_12m")
### discrete outcomes
### zeros correspond to utilization and heath (trimming from below)
### ones correspond to cost outcomes (trimming from above)

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
                                       stratification_controls)


### weight
weight_name<-"weight_12m"  
### quantile grid size for quantile regression
quantile_grid_size<-0.01

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

formula_threshold_list[["numeric_health"]]<-as.formula(paste0("selection~treat+( ",paste0(list_of_all_controls,collapse="+"),")" ))
formula_point_mass_list[["numeric_health"]]<-as.formula(paste0("outcome_point_mass~ ",paste0(stratification_controls,collapse="+") ))

formula_threshold_list[["cost"]]<-as.formula(paste0("selection~treat*(female_list
                                                     +zip_msa+snap_tot_hh_prenotify07+
                                                     ed_charg_tot_pre_ed+snap_ever_prenotify07+
                                                     charg_tot_pre_ed+num_epct_pre_ed+
                                                     num_out_pre_cens_ed+num_visit_pre_cens_ed+
                                                     any_out_pre_ed+num_edcnnp_pre_ed+
                                                     any_hiun_pre_ed+num_hiun_pre_cens_ed+
                                                     num_ne_pre_ed+num_on_pre_cens_ed+
                                                     tanf_tot_hh_prenotify07+any_on_pre_ed+
                                                     any_visit_pre_ed)*(",paste0(list_of_all_controls,collapse="+"),")" ))



prepared_data_sample_12m<-sapply(prepared_data_sample_12m,as.numeric)
prepared_data_sample_12m<-as.data.frame(prepared_data_sample_12m)

trimmed_data<-list()
result<-list()
n_trimmed<-list()





### continuous outcomes with point mass at trimming value
for (outcome_name in c(utilization_outcomes)) {
  print(outcome_name)
  trimmed_data[[outcome_name]]<-first_stage_wrapper_continuous(data=prepared_data_sample_12m,
                                                             outcome_name=outcome_name,
                                                             list_of_all_controls=list_of_all_controls,
                                                             selection_function_for_threshold=glm,
                                                             selection_function_for_threshold_name="glm",
                                                             stratification_controls= stratification_controls,
                                                             selection_function_point_mass=rlassologit,
                                                             selection_function_point_mass_name="rlassologit",
                                                             formula_threshold=formula_threshold_list[["utilization"]],
                                                             formula_point_mass=formula_point_mass_list[["utilization"]],
                                                             myweights=prepared_data_sample_12m$weight_12m,
                                                             variables_for_outcome=selected_covariates_all_regressions,
                                                             quantile_grid_size=quantile_grid_size,
                                                             sort_quantiles=TRUE,
                                                             selected_inds=NULL,
                                                             max_p0_hat=1.001001,
                                                             point_mass_cutoff=0,
                                                             point_mass=FALSE,
                                                             trimming_value=0)
  n_trimmed[[outcome_name]]<-dim(trimmed_data[[outcome_name]])[1]
  result[[outcome_name]]<-second_stage_wrapper(data=trimmed_data[[outcome_name]],
                                               weight_name=weight_name,
                                               stratification_controls=stratification_controls)
}

for (outcome_name in c(numeric_health_outcomes)) {
  trimmed_data[[outcome_name]]<-first_stage_wrapper_continuous(data=prepared_data_sample_12m,
                                                                          outcome_name=outcome_name,
                                                                          list_of_all_controls=list_of_all_controls,
                                                                          selection_function_for_threshold=glm,
                                                                          selection_function_for_threshold_name="glm",
                                                                          stratification_controls= stratification_controls,
                                                                          selection_function_point_mass=rlassologit,
                                                                          selection_function_point_mass_name="rlassologit",
                                                                          formula_threshold=formula_threshold_list[["numeric_health"]],
                                                                          formula_point_mass=formula_point_mass_list[["numeric_health"]],
                                                                          myweights=prepared_data_sample_12m$weight_12m,
                                                                          variables_for_outcome=stratification_controls,
                                                                          quantile_grid_size=0.01,
                                                                          sort_quantiles=TRUE,
                                                                          selected_inds=NULL,
                                                                          max_p0_hat=1.001001,
                                                                          point_mass_cutoff=0.18,
                                                               point_mass=TRUE, trimming_value=0)
  print(dim( trimmed_data[[outcome_name]]))
  result[[outcome_name]]<-second_stage_wrapper(data=trimmed_data[[outcome_name]],
                                               weight_name=weight_name,
                                               stratification_controls=stratification_controls)
}

for (outcome_name in discrete_outcome_names) {
  
  
  trimmed_data[[outcome_name]]<-first_stage_wrapper_discrete(data=prepared_data_sample_12m,
                                                             outcome_name=outcome_name,
                                                             list_of_all_controls=list_of_all_controls,
                                                             trimming_value=trimming_value_list[match(outcome_name,discrete_outcome_names)],
                                                             selection_function_for_threshold=rlassologit,
                                                             selection_function_for_threshold_name="rlassologit",
                                                             stratification_controls= stratification_controls,
                                                             selection_function_point_mass=rlassologit,
                                                             selection_function_point_mass_name="rlassologit",
                                                             formula_threshold=formula_threshold_list[["discrete"]],
                                                             formula_point_mass=formula_point_mass_list[["discrete"]],
                                                             myweights=prepared_data_sample_12m$weight_12m,
                                                             point_mass_cutoff=0.10,
                                                             max_p0_hat=1.000001 )
  
  
  
  
  result[[outcome_name]]<-second_stage_wrapper(data=trimmed_data[[outcome_name]],
                                               weight_name=weight_name,
                                               stratification_controls=stratification_controls)
  
}


trim_ml<-matrix(unlist(result),nrow=3,byrow=FALSE)
rownames(trim_ml)<-c("control_mean","itt","late")
colnames(trim_ml)<-setdiff(all_names,cost_outcomes)
trim_ml<-t(trim_ml)

#basic_trim_sd<-basic_trim_sd[,match(colnames(basic_trim_sd),all_names)]
#colnames(basic_trim_sd)<-all_names

write.csv(trim_ml,paste0(my_path,"/OHIE/STEP3_ML/csv/estimates_ml_trim.csv"))
#write.csv(basic_trim_sd,paste0(my_path,"/OHIE/STEP3_ML/csv/sd_ml_trim.csv"))


save.image(paste0(my_path,"/OHIE/STEP3_ML/draft/estimates.RData"))
