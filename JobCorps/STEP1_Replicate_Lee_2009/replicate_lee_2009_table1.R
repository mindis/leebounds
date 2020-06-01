rm(list=ls())

library(tidyverse)
library(radiant.data)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
my_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
my_data<-as.data.frame(my_data)

varnames<-c("FEMALE","AGE","WHITE","BLACK","HISP","OTHERRAC","NEVERMARRIED",
            "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
            "HH_INC1","HH_INC2","HH_INC3","HH_INC4","HH_INC5",
            "PERS_INC1","PERS_INC2","PERS_INC3","PERS_INC4",
            
            "CURRJOB.x", "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
            "WKEARNR",
            "EARNQ4","EARNQ8","EARNQ12","EARNQ16",
            
            "HWH52","HWH104","HWH156","HWH208","EARNH52","EARNH104","EARNH156","EARNH208"
)

my_data_treated<-filter(my_data,TREATMNT.y==1)
my_data_control<-filter(my_data,TREATMNT.y==0)

dim(my_data_control)[1]

dim(my_data_treated)[2]
control_prop_pres<-sapply(my_data_control[,varnames],function (x) weighted.mean(!is.na(x),w=my_data_control$DSGN_WGT.y))
control_mean<-sapply(my_data_control[,varnames],weighted.mean,w=my_data_control$DSGN_WGT.y,na.rm=TRUE)
control_sd<-sapply(my_data_control[,varnames],weighted.sd,w=my_data_control$DSGN_WGT.y,na.rm=TRUE)

treated_prop_pres<-sapply(my_data_treated[,varnames],function (x) weighted.mean(!is.na(x),w=my_data_treated$DSGN_WGT.y))
treated_mean<-sapply(my_data_treated[,varnames],weighted.mean,w=my_data_treated$DSGN_WGT.y,na.rm=TRUE)
treated_sd<-sapply(my_data_treated[,varnames],weighted.sd,w=my_data_treated$DSGN_WGT.y,na.rm=TRUE)

table1<-cbind(control_prop_pres,control_mean,control_sd,
              treated_prop_pres,treated_mean,treated_sd)
colnames(table1)<-c("CONTROL_PRESENT","CONTROL_AVE","CONTROL_SD","PROGRAM_PRESENT","PROGRAM_AVE","PROGRAM_SD")

table1<-apply(table1,2,round,3)
rownames(table1)<-varnames

write.csv(table1,paste0(my_path,"/JobCorps/STEP1_Replicate_Lee_2009/csv/table1.csv"))
