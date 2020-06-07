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

## for parallel computing
library(doParallel)
library(foreach)
cl <- makeCluster(8)
registerDoParallel(cl)
on.exit(stopCluster(cl))

my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path, "/OHIE/STEP2_Standard_Trim/utils.R"))
########### READ IN DATA ################

prepared_data<-read_feather(paste0(my_path,"/OHIE/OHIE_Data/prepared_data.feather"))
prepared_data_sample_12m<-filter(prepared_data,sample_12m=="In 12m mail survey sample")
prepared_data_responded_control<-filter(prepared_data_sample_12m,treatment == 0)

### Unique covariate groups determined by compulsory covariates 
print(round(table(prepared_data_sample_12m$group_id)/dim(prepared_data_sample_12m)[1]*100,3))
### Merge groups 16,17,18 into one
prepared_data_sample_12m$group_id [prepared_data_sample_12m$group_id %in% c(16,17,18)]<-16
unique_groups<-unique(prepared_data_sample_12m$group_id)
prepared_data_sample_12m$group_id[prepared_data_sample_12m$group_id==16]<-3
stratification_controls<-grep("ddd",colnames(prepared_data),value=TRUE)
weight_name<-"weight_12m"
############# OUTCOMES #########


discrete_outcome_names<-c("rx_any_12m","doc_any_12m","er_any_12m","hosp_any_12m",
                          "cost_any_oop_12m","cost_any_owe_12m","cost_borrow_12m","cost_refused_12m", 
                          "health_genflip_bin_12m","health_notpoor_12m","health_chgflip_bin_12m","nodep_screen_12m")

utilization_outcomes<-c("rx_num_mod_12m","doc_num_mod_12m","er_num_mod_12m","hosp_num_mod_12m"
)

cost_outcomes<-c(   "cost_doc_oop_12m",        "cost_er_oop_12m"  ,       "cost_rx_oop_12m",         "cost_oth_oop_12m" ,       "cost_tot_oop_12m" ,      
                    "cost_tot_owe_12m" ,          "cost_tot_oop_mod_12m" ,   "cost_tot_owe_mod_12m")

numeric_health_outcomes<-c("notbaddays_phys_12m","notbaddays_tot_12m","notbaddays_ment_12m")
all_names<-c(utilization_outcomes,cost_outcomes,numeric_health_outcomes,discrete_outcome_names)




prepared_data_sample_12m<-prepared_data_sample_12m[,c("treatment",stratification_controls,
                                                      "group_id","ohp_all_ever_survey","household_id","weight_12m",
                                                      all_names,"person_id")]

weight_name<-"weight_12m"
B<-1000
############# Estimate bounds ############

### trimming value is zero for utilization or health outcomes
trimming_value_all_names<-rep(0,length(all_names))
## trimming value is 1 for discrete cost outcomes
trimming_value_all_names[grep("cost",all_names)]<-1

flag_discrete_all_names<-rep(TRUE,length(all_names))
flag_discrete_all_names[13:15]<-FALSE
flag_discrete_all_names[5:12]<-FALSE


cluster_robust_se<-list()
control_sample_size<-list()

table_basic_trim<-read.csv(paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_basic_trim.csv"),row.names=1)
table_no_trim<-read.csv(paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_no_trim.csv"),row.names=1)

for (i in 1:length(all_names)) {
  outcome_name =all_names[i]
  print(outcome_name)
  myres<- foreach(b=1:B, .combine = 'cbind',.packages = c("tidyverse","AER","miceadds"))  %dopar% {
    ## some (a small number) bootstrap draws may fail the basic trimming method
    ans<-tryCatch(bootstrap_mean_itt_late(b,outcome_name=outcome_name,weight_name="weight_12m",
                                          trimming_value=trimming_value_all_names[match(outcome_name,all_names)],
                                          flag_discrete=flag_discrete_all_names[match(outcome_name,all_names)]),error=function(e) NULL)
    
  }
  
  
  cluster_robust_se[[outcome_name]]<-apply(myres,1,sd)
 # control_sample_size[[outcome_name]]<-sqrt(sum(!is.na(prepared_data_responded_control[,outcome_name])))
 
  
 # cluster_robust_se[[outcome_name]][1]<- cluster_robust_se[[outcome_name]][1]*control_sample_size[[outcome_name]]
 
}
basic_trim_sd<-matrix(unlist(cluster_robust_se),nrow=4,byrow=FALSE)
rownames(basic_trim_sd)<-c("control_mean","itt","late","control_mean_notrim")
#colnames(basic_trim_sd)<-c(utilization_outcomes,numeric_health_outcomes,discrete_outcome_names[1:4],discrete_outcome_names[9:12])
colnames(basic_trim_sd)<-all_names
basic_trim_sd<-t(basic_trim_sd)
basic_trim_sd<-apply(basic_trim_sd,2,round,4)

write.csv(basic_trim_sd,paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/sd_basic_trim_clustered.csv"))




#save.image(paste0(my_path,"/OHIE/STEP2_Standard_Trim/draft/std_trim_boot.RData"))

