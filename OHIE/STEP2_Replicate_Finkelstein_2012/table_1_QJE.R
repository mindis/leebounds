#### This code replicates selected rows of Table 1, Finkelstein et al (2014), QJE, page 1068-1069
###   by Vira Semenova

rm(list=ls())
#install.packages("sandwich",lib="~/apps/R_3.5.1/")
library(foreign,lib.loc="~/apps/R_3.5.1/")
library(tidyverse,lib.loc="~/apps/R_3.5.1/")



########### READ IN DATA ################
my_path<-"/net/holyparkesec/data/tata/leebounds/"
baseline<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/oregonhie_descriptive_vars.dta"))
survey12<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/oregonhie_survey12m_vars.dta"))
prepared_data<-read.dta(paste0(my_path,"/OHIE/OHIE_Data/data_for_analysis_ver12.dta"))

survey12<-left_join(survey12, baseline[, c("treatment","person_id","female_list","english_list","zip_msa_list")],by="person_id") 
survey12<-left_join(survey12, prepared_data)
baseline<-left_join(baseline,prepared_data)
survey12_responded<-survey12[survey12$sample_12m_resp == "12m mail survey responder",]


########## REPLICATE TABLE 1 (MAIN TEXT QJE) PANEL A ###############
baseline_control<-filter(baseline,treatment == "Not selected")


##  ccontrol mean sex
baseline_control$female_list<-as.character(baseline_control$female_list)
mean(baseline_control$female_list=="1: Female")
##   control mean age 50-64
mean(baseline_control$older,na.rm=TRUE)
##  control mean age 20-50
mean(baseline_control$younger,na.rm=TRUE)
##  control mean english preferred language
mean(baseline_control$english_list == "Requested English materials")
##  zip code in MSA
mean(baseline_control$zip_msa_list == "Zip code of residence in a MSA",na.rm=TRUE)




################################# PANEL B ###############
survey12_responded_control<-filter(survey12_responded,treatment == "Not selected" & sample_12m_resp=="12m mail survey responder")

weighted.mean(survey12_responded_control$female_list =="1: Female",w=survey12_responded_control$weight_12m,na.rm=TRUE)

weighted.mean(survey12_responded_control$older,na.rm=TRUE,w=survey12_responded_control$weight_12m)
##  control mean age 20-50
weighted.mean(survey12_responded_control$younger,na.rm=TRUE,w=survey12_responded_control$weight_12m)

weighted.mean(survey12_responded_control$english_list == "Requested English materials" ,na.rm=TRUE,w=survey12_responded_control$weight_12m)

weighted.mean(survey12_responded_control$zip_msa_list == "Zip code of residence in a MSA" ,na.rm=TRUE,w=survey12_responded_control$weight_12m)

############# PANEL B: EDUCATION

weighted.mean(survey12_responded_control$edu_12m=="less than hs",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$edu_12m=="vocational or 2-year degree",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$edu_12m=="4-year degree",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$edu_12m=="hs diploma or GED",na.rm=TRUE,w=survey12_responded_control$weight_12m)

############# PANEL B: EMPLOYMENT

weighted.mean(survey12_responded_control$employ_hrs_12m=="don't currently work",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$employ_hrs_12m=="work <20 hrs/week",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$employ_hrs_12m=="work 20-29 hrs/week",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$employ_hrs_12m=="work 30+ hrs/week",na.rm=TRUE,w=survey12_responded_control$weight_12m)

############# PANEL B: HOUSEHOLD INCOME

weighted.mean(survey12_responded_control$hhinc_pctfpl_12m<=50,na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$hhinc_pctfpl_12m>=50 &
                survey12_responded_control$hhinc_pctfpl_12m<75,na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$hhinc_pctfpl_12m>=75 &
                survey12_responded_control$hhinc_pctfpl_12m<100,na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$hhinc_pctfpl_12m>=100 &
                survey12_responded_control$hhinc_pctfpl_12m<150,na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$hhinc_pctfpl_12m>=150,na.rm=TRUE,w=survey12_responded_control$weight_12m)

############# PANEL B: INSURANCE

weighted.mean(survey12_responded_control$ins_any_12m=="Yes",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$ins_ohp_12m=="Yes",na.rm=TRUE,w=survey12_responded_control$weight_12m)
weighted.mean(survey12_responded_control$ins_private_12m=="Yes",na.rm=TRUE,w=survey12_responded_control$weight_12m)



