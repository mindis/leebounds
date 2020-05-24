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
list_of_basic_qje_controls<-grep("ddd",colnames(prepared_data),value=TRUE)
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

### INTERMEDIATE AFTER TRIMMING SUMMARY STAT

trimmed_data<-list()

### final treat-control difference after trimming
delta_after_trimming<-matrix(0,length(unique(prepared_data_sample_12m$group_id)),length(all_names))
n_trimmed<-matrix(0,length(unique_groups),length(all_names))
prop_trimmed<-matrix(0,length(unique_groups),length(all_names))
prop_trimmed_lowest_value<-matrix(0,length(unique_groups),length(all_names))



############ SAVE FINAL RESULTS HERE ############
table_control_means_notrim<-matrix(0,length(all_names),2)
rownames(table_control_means_notrim)<-all_names
colnames(table_control_means_notrim)<-c("mean","sd")

table_itt_notrim<-matrix(0,length(all_names),2)
rownames(table_itt_notrim)<-all_names
colnames(table_itt_notrim)<-c("mean","sd")

table_late_notrim<-matrix(0,length(all_names),2)
rownames(table_late_notrim)<-all_names
colnames(table_late_notrim)<-c("mean","sd")


table_control_means_basic_trim<-matrix(0,length(all_names),2)
rownames(table_control_means_basic_trim)<-all_names
colnames(table_control_means_basic_trim)<-c("mean","sd")

table_itt_basic_trim<-matrix(0,length(all_names),2)
rownames(table_itt_basic_trim)<-all_names
colnames(table_itt_basic_trim)<-c("mean","sd")

table_late_basic_trim<-matrix(0,length(all_names),2)
rownames(table_late_basic_trim)<-all_names
colnames(table_late_basic_trim)<-c("mean","sd")



############# Estimate bounds ############

### trimming value is zero for utilization or health outcomes
trimming_value_all_names<-rep(0,length(all_names))
## trimming value is 1 for discrete cost outcomes
trimming_value_all_names[grep("cost",all_names)]<-1

flag_discrete_all_names<-rep(TRUE,length(all_names))
flag_discrete_all_names[13:15]<-FALSE
flag_discrete_all_names[5:12]<-FALSE

for (i in 1:length(all_names)) {
  outcome_name =all_names[i]
  trimmed_data[[i]]<-matrix(0,nrow=0,ncol=dim(prepared_data))
  
  for (j in 1:length(unique(prepared_data_sample_12m$group_id)) ) {
    res<-basic_trimming_subjects(data=prepared_data_sample_12m,
                                 outcome_name= outcome_name, 
                                 group_ind=j,
                                 trimming_value=trimming_value_all_names[i],
                                 step=1,
                                 flag_discrete = flag_discrete_all_names[i])
    
    control_sample=res$control_sample
    treated_sample=res$treated_sample
    nontrimmed_controls<-res$nontrimmed_controls
    trimmed_data[[i]]<-rbind(trimmed_data[[i]],
                             control_sample[ control_sample$person_id %in% nontrimmed_controls,],
                             treated_sample)
    delta_after_trimming[j,i]<-res$delta_after_trimming
    n_trimmed[j,i]<-res$n_trimmed
    prop_trimmed[j,i]<-res$prop_trimmed
    prop_trimmed_lowest_value[j,i]<-res$prop_trimmed_lowest_value
    
  } 
  
}

print(all(delta_after_trimming>=0))

##### check that response rate of the controls and the treated are approximately equal
### compute maximum deviation
print(max(abs(delta_after_trimming)))


prepared_data_sample_12m<-as.data.frame(prepared_data_sample_12m)
prepared_data_responded_control<-as.data.frame(prepared_data_responded_control)
for (i in 1:length(all_names)) {
  outcome_name =all_names[i]
  ## column 1 
  table_control_means_notrim[i,1]<-round(weighted.mean(prepared_data_responded_control[,outcome_name] ,na.rm=TRUE,w=prepared_data_responded_control[,weight_name]),3)
  table_control_means_notrim[i,2]<-round(radiant.data::weighted.sd(prepared_data_responded_control[,outcome_name] ,na.rm=TRUE,w=prepared_data_responded_control[,weight_name]),3)
  
  ## column 5 
  control_sample<-as.data.frame(trimmed_data[[i]][trimmed_data[[i]]$treatment==0,])
  table_control_means_basic_trim[i,1]<-round(weighted.mean(control_sample[,outcome_name] ,na.rm=TRUE,w=control_sample[, weight_name]),3)
  table_control_means_basic_trim[i,2]<-round(radiant.data::weighted.sd(control_sample[,outcome_name] ,na.rm=TRUE,w=control_sample[, weight_name] ),3)
  
}
table_control_means_notrim
table_control_means_basic_trim

################# REPLICATE COLUMN 2 and 5 (ITT) ###################


for (i in 1:length(all_names)) {
  outcome_name =all_names[i]
  form<-paste0(outcome_name,"~treatment+",paste0(list_of_basic_qje_controls,collapse="+"))
  
  ## column 2
  res<-miceadds::lm.cluster(formula=as.formula(form), data=prepared_data_sample_12m,weights=prepared_data_sample_12m[,weight_name],
                            cluster="household_id")
  table_itt_notrim[i,1]<-round(res$lm_res$coefficients[2],4)
  table_itt_notrim[i,2]<-round(summary(res)[2,2],4)
  
  ## column 6
  trimmed_data[[i]]<-as.data.frame(trimmed_data[[i]])
  res<-miceadds::lm.cluster(formula=as.formula(form), data=trimmed_data[[i]],weights=trimmed_data[[i]][,weight_name],
                            cluster="household_id")
  table_itt_basic_trim[i,1]<-round(res$lm_res$coefficients[2],4)
  table_itt_basic_trim[i,2]<-round(summary(res)[2,2],4)
}

################# REPLICATE COLUMN 3 and 6 (LATE) ###################


for (i in 1:length(all_names)) {
  outcome_name =all_names[i]
  ## column 3
  form<-paste0( paste0(outcome_name,"~ohp_all_ever_survey+",paste0(list_of_basic_qje_controls,collapse="+")),"|treatment+",paste0(list_of_basic_qje_controls,collapse="+"))
  res<-AER::ivreg(formula=as.formula(form), data=prepared_data_sample_12m[!is.na(prepared_data_sample_12m[,outcome_name]),],
                  weights=prepared_data_sample_12m[!is.na(prepared_data_sample_12m[,outcome_name]),weight_name]  )
  
  table_late_notrim[i,1]<-round(res$coefficients[2],4)
  cluster.robust<-ivpack::cluster.robust.se(res,clusterid=prepared_data_sample_12m$household_id[!is.na(prepared_data_sample_12m[,outcome_name])])
  table_late_notrim[i,2]<-round(cluster.robust[2,2],4)
  
  ## column 7
  trimmed_data[[i]]<-as.data.frame(trimmed_data[[i]])
  res<-AER::ivreg(formula=as.formula(form), data=trimmed_data[[i]][!is.na(trimmed_data[[i]][,outcome_name]),],weights=trimmed_data[[i]][!is.na(trimmed_data[[i]][,outcome_name]),weight_name]  )
  
  table_late_basic_trim[i,1]<-round(res$coefficients[2],4)
  cluster.robust<-ivpack::cluster.robust.se(res,clusterid=trimmed_data[[i]]$household_id[!is.na(trimmed_data[[i]][,outcome_name])])
  table_late_basic_trim[i,2]<-round(cluster.robust[2,2],4)
}

### save results
table_estimates_notrim<-data.frame(control_mean=table_control_means_notrim[,1],
                                   itt_mean=table_itt_notrim[,1],
                                   late_mean=table_late_notrim[,1])

table_sd_notrim<-data.frame(control_mean_sd=table_control_means_notrim[,2],
                            itt_sd=table_itt_notrim[,2],
                            late_sd=table_late_notrim[,2])

table_notrim<-cbind(table_estimates_notrim,table_sd_notrim)
table_estimates_basic_trim<-data.frame(control_mean=table_control_means_basic_trim[,1],
                                       itt_mean=table_itt_basic_trim[,1],
                                       late_mean=table_late_basic_trim[,1])
#non-trimming robust standard error
## computed for reference 
## correct (trimming-robust) standard error is computed in bootstrap_lower_bound_std_trim.R
table_sd_basic_trim<-data.frame(control_mean_sd=table_control_means_basic_trim[,2],
                                itt_sd=table_itt_basic_trim[,2],
                                late_sd=table_late_basic_trim[,2])
table_basic_trim<-cbind(table_estimates_basic_trim,table_sd_basic_trim)

write.csv(table_notrim,paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_no_trim.csv"))
write.csv(table_basic_trim,paste0(my_path,"/OHIE/STEP2_Standard_Trim/csv/estimates_basic_trim.csv"))


