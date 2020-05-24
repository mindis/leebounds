# clear workspace
rm(list=ls())

library(foreign)
library(tidyverse)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
prepared_data<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/data_for_analysis_ver12.dta"))
ed_vars<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/oregonhie_ed_vars.dta"))
state_vars<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/oregonhie_stateprograms_vars.dta"))


## this is treatment variable (LOTTERY)
prepared_data$treatment<-as.character(prepared_data$treatment)
prepared_data$treatment[prepared_data$treatment=="Not selected"]<-0
prepared_data$treatment[prepared_data$treatment=="Selected"]<-1
prepared_data$treatment<-as.numeric(prepared_data$treatment)

## this is instrument variable (INSURANCE)
prepared_data$ohp_all_ever_survey<-as.character(prepared_data$ohp_all_ever_survey)
prepared_data$ohp_all_ever_survey[prepared_data$ohp_all_ever_survey=="Enrolled"]<-"1"
prepared_data$ohp_all_ever_survey[prepared_data$ohp_all_ever_survey=="NOT enrolled"]<-"0"
prepared_data$ohp_all_ever_survey<-as.numeric(prepared_data$ohp_all_ever_survey)

## these are compulsory covariates  X (used for stratification)
grep("ddd",colnames(prepared_data),value=TRUE)

## these are other baseline covariates (all covariates that have _pre_ in their name)
ed_pre_covs<-grep("_pre_",colnames(ed_vars),value=TRUE)
state_pre_covs<-grep("prenotify",colnames(state_vars),value=TRUE)

prepared_data<-left_join(prepared_data,ed_vars[,c(ed_pre_covs,"person_id")])
prepared_data<-left_join(prepared_data,state_vars[,c(state_pre_covs,"person_id")])

################## MAIN OUTCOME VARIABLES OF INTEREST #######


### binary outcomes (extensive margin)
binary_outcome_names<-c("rx_any_12m","doc_any_12m","er_any_12m","hosp_any_12m","cost_any_oop_12m","cost_any_owe_12m","cost_borrow_12m","cost_refused_12m",
                        "health_genflip_bin_12m","health_notpoor_12m","health_chgflip_bin_12m","nodep_screen_12m"  )


#### continuous outcomes (total utilization)
utilization_outcomes<-c("rx_num_mod_12m","doc_num_mod_12m","er_num_mod_12m","hosp_num_mod_12m")

#### cost outcomes
cost_outcomes<-c(   "cost_doc_oop_12m",        "cost_er_oop_12m"  ,       "cost_rx_oop_12m",         "cost_oth_oop_12m" ,       "cost_tot_oop_12m" ,      
                 "cost_tot_owe_12m" ,                  "cost_tot_oop_mod_12m" ,   "cost_tot_owe_mod_12m")

numeric_health_outcomes<-c("notbaddays_phys_12m","notbaddays_tot_12m","notbaddays_ment_12m")
all_names<-c(binary_outcome_names,utilization_outcomes,cost_outcomes,numeric_health_outcomes)
################## STEP 1: BASIC CLEANING #######################



prepared_data$age<-2008-prepared_data$birthyear_list

replace_yes_no<-function(x) {
  x<-as.character(x)
  x[x=="Yes"]<-"1"
  x[x=="No"]<-"0"
  return(x)
}


#prepared_data<-sapply(prepared_data,as.character)   

prepared_data[,grep("any_|_ever_",c(ed_pre_covs,state_pre_covs),value=TRUE)]<-sapply(prepared_data[,grep("any_|_ever_",c(ed_pre_covs,state_pre_covs),value=TRUE)],
                                                                                     replace_yes_no)

prepared_data$female_list<-as.character(prepared_data$female_list)
prepared_data$female_list[prepared_data$female_list=="1: Female"]<-"1"
prepared_data$female_list[prepared_data$female_list=="0: Male"]<-"0"

prepared_data$english_list<-as.character(prepared_data$english_list)
prepared_data$english_list[prepared_data$english_list=="Requested English materials"]<-"1"
prepared_data$english_list[prepared_data$english_list=="Requested materials in language other than english"]<-"0"

prepared_data$zip_msa<-as.character(prepared_data$zip_msa)
prepared_data$zip_msa[prepared_data$zip_msa=="Zip code of residence in a MSA"]<-"1"
prepared_data$zip_msa[prepared_data$zip_msa=="Zip code of residence NOT in a MSA"]<-"0"



prepared_data$age<-2008-prepared_data$birthyear_list




prepared_data[,binary_outcome_names]<-sapply(prepared_data[,binary_outcome_names],as.character)
prepared_data[,binary_outcome_names]<-sapply(prepared_data[,binary_outcome_names],replace_yes_no)
prepared_data[,all_names]<-sapply(prepared_data[,all_names],as.numeric)



write.csv(prepared_data[,c("person_id","household_id","treatment", "ohp_all_ever_survey","female_list","english_list","zip_msa",
                               ed_pre_covs,state_pre_covs,
                               grep("ddd",colnames(prepared_data),value=TRUE),
                               "weight_12m","sample_12m","sample_12m_resp",
                               binary_outcome_names,
                               utilization_outcomes,
                               cost_outcomes,numeric_health_outcomes
)],paste0(my_path,"/OHIE/OHIE_Data/prepared_data.csv"))

write_feather(prepared_data[,c("person_id","household_id","treatment", "ohp_all_ever_survey","female_list","english_list","zip_msa",
                               ed_pre_covs,state_pre_covs,
                               grep("ddd",colnames(prepared_data),value=TRUE),
                               "weight_12m","sample_12m","sample_12m_resp",
                               binary_outcome_names,
                               utilization_outcomes,
                               cost_outcomes,numeric_health_outcomes
)],paste0(my_path,"/OHIE/OHIE_Data/prepared_data.feather"))
