rm(list=ls())
## load libraries

library(sas7bdat)
library(tidyverse)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
fulldata_2002<-read.sas7bdat(paste0(my_path,"/ABBKK_2002/data/aerdat4.sas7bdat"))
tab5_test_scores<-read.sas7bdat(paste0(my_path,"/ABBKK_2002/data/tab5v1.sas7bdat"))
fulldata_2002_scores<-left_join(fulldata_2002,tab5_test_scores[,c("ID","DAD_MISS", "MOM_MISS", "MATH",     "READING",  "WRITING","TOTALPTS",
                                                                  "TSITE1",   "TSITE2" ,  "TSITE3"  )],
                                by=c("ID"="ID"))
#### HANDLING MISSING VALUES IN COVARIATES #################
fulldata_2002_scores<-filter(fulldata_2002_scores,BOG95==1  & (!is.na(VOUCH0)) & !is.na(AGE2))
fulldata_2002_scores$MOM_AGE_IS_NA<-is.na(fulldata_2002_scores$MOM_AGE)
fulldata_2002_scores$DAD_AGE_IS_NA<-is.na(fulldata_2002_scores$DAD_AGE)
fulldata_2002_scores$MOM_SCH_IS_NA<-is.na(fulldata_2002_scores$MOM_SCH)
fulldata_2002_scores$DAD_SCH_IS_NA<-is.na(fulldata_2002_scores$DAD_SCH)
fulldata_2002_scores$MOM_MW_IS_NA<-is.na(fulldata_2002_scores$MOM_MW)
fulldata_2002_scores$DAD_MW_IS_NA<-is.na(fulldata_2002_scores$DAD_MW)
fulldata_2002_scores$SEX_IS_NA<-is.na(fulldata_2002_scores$SEX_NAME)



fulldata_2002_scores$MOM_AGE[is.na(fulldata_2002_scores$MOM_AGE)]<-0
fulldata_2002_scores$DAD_AGE[is.na(fulldata_2002_scores$DAD_AGE)]<-0
fulldata_2002_scores$MOM_SCH[is.na(fulldata_2002_scores$MOM_SCH)]<-0
fulldata_2002_scores$DAD_SCH[is.na(fulldata_2002_scores$DAD_SCH)]<-0
fulldata_2002_scores$MOM_MW[is.na(fulldata_2002_scores$MOM_MW)]<-0
fulldata_2002_scores$DAD_MW[is.na(fulldata_2002_scores$DAD_MW)]<-0
fulldata_2002_scores$SEX_NAME[is.na(fulldata_2002_scores$SEX_NAME)]<-0

fulldata_2002_scores$MATH[is.na(fulldata_2002_scores$MATH)]<-0
fulldata_2002_scores$READING[is.na(fulldata_2002_scores$READING)]<-0
fulldata_2002_scores$WRITING[is.na(fulldata_2002_scores$WRITING)]<-0

### dropped DAREA and STRATA codes that had <0.001 non-zero cases in data
## example :DAREA3
summary(fulldata_2002_scores$DAREA3)
exogenous_covariates<-c("AGE2","SEX_NAME","MOM_SCH" ,"MOM_AGE","DAD_SCH" ,"DAD_AGE",
                        "DAREA4",   "DAREA5",   "DAREA6",   "DAREA7",
                        "DAREA11" ,  "DAREA15",
                        "DAREA16" , "DAREA17"  ,"DAREA18" , "DAREA19","MOM_SCH_IS_NA","MOM_AGE_IS_NA","DAD_SCH_IS_NA","DAD_AGE_IS_NA",
                        "SEX_IS_NA",
                        "STRATA1",  "STRATA2",  "STRATA3" , "STRATA4","STRATA5","STRATA6",
                        "DAD_MW","MOM_MW","DAD_MW_IS_NA","MOM_MW_IS_NA",
                        "TSITE1",   "TSITE2" ,  "TSITE3"  )
# wages are collected at survey conducted after random assignment / may be correlated with voucher


## 3. Define the data set
mydata<-fulldata_2002_scores[,c("VOUCH0","TEST_TAK","MATH","READING","WRITING",exogenous_covariates)]
sum(is.na(mydata))
# same prop of vouch0 as in ABBKK (2002)
mean(mydata$VOUCH0)
mydata[mydata==TRUE]<-1
mydata[mydata==FALSE]<-0
write.csv(mydata,paste0(my_path,"/ABBKK_2002/data/ABK_prepared_data.csv"))
