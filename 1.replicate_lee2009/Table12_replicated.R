## read in data
rm(list=ls())
library(sas7bdat)
library(dplyr)
library(tidyverse)
library(SDMTools)
library(feather)
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"~"
my_data<-read_feather(paste0(my_path,"/R_code/dataLee2009.feather"))
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


average_cols_weighted<-matrix(0,length(varnames),2)
sd_cols_weighted<-matrix(0,length(varnames),2)
prop_na_weighted<-matrix(0,length(varnames),2)

weight_0<-my_data$DSGN_WGT.y[my_data$TREATMNT.y==0]
  
weight_0<-  weight_0/sum(weight_0)
weight_1<-my_data$DSGN_WGT.y[my_data$TREATMNT.y==1]
weight_1<-  weight_1/mean(weight_1)
weights<-list(weight_0,weight_1)

for (i in 1:length(varnames)) {
  for (treat in c(0,1)) {
    varname<-varnames[i]
    weight_treat<-weights[[treat+1]]
    average_cols_weighted[i,treat+1]=wt.mean(x=my_data[my_data$TREATMNT.y==treat,varname],
                                             wt=weight_treat)
    sd_cols_weighted[i,treat+1]=wt.sd(x=my_data[my_data$TREATMNT.y==treat,varname],
                                      wt=weight_treat)
    prop_na_weighted[i,treat+1]=wt.mean(!is.na(my_data[my_data$TREATMNT.y==treat,varname]),
                                        wt=weight_treat)
  }
}
rownames(average_cols_weighted)<-varnames
rownames(sd_cols_weighted)<-varnames
rownames(prop_na_weighted)<-varnames
colnames(average_cols_weighted)<-c("CONTROL","PROGRAM")
colnames(prop_na_weighted)<-c("CONTROL","PROGRAM")
colnames(sd_cols_weighted)<-c("CONTROL","PROGRAM")

average_cols_weighted<-round(average_cols_weighted,3)
sd_cols_weighted<-round(sd_cols_weighted,3)
prop_na_weighted<-round(prop_na_weighted,3)

table1<-cbind(average_cols_weighted,sd_cols_weighted,prop_na_weighted)
colnames(table1)<-c("CONTROL_AVE","PROGRAM_AVE","CONTROL_SD","PROGRAM_SD","CONTROL_PRESENT","PROGRAM_PRESENT")
write.csv(table1,paste0(my_path,"/R_code/1.replicate_lee2009/tables/Table1_replicated.csv"))


### Replication of Table 2 Lee (2009)


my_data$EMPLOYED_208<-as.numeric(my_data$HWH208>0)
my_data$EMPLOYED_208<-as.factor(my_data$EMPLOYED_208)
int_weight<-round(my_data$DSGN_WGT.y,0)
glm_model_weighted<-glm(EMPLOYED_208~.,data=my_data[,c("TREATMNT.y","FEMALE","AGE","BLACK","HISP","OTHERRAC",
                                                              "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                                                              "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                                                          "PERS_INC2","PERS_INC3","PERS_INC4",
                                                              
                                                              "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                                                              "WKEARNR","EMPLOYED_208"
)],weights=int_weight,family = binomial)



coefs<-round(coefficients(glm_model_weighted),6)
se<-round(coefficients(glm_model_weighted),6)
table2=round(summary(glm_model_weighted)$coefficients[,1:2],3)
write.csv(table2,paste0(my_path,"/R_code/1.replicate_lee2009/tables/Table2_replicated.csv"))
