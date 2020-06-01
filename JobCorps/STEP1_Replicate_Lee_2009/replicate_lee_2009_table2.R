## read in data
rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")
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

my_data$EMPLOYED_208<-as.numeric(my_data$HWH208>0)
my_data$EMPLOYED_208<-as.factor(my_data$EMPLOYED_208)


#int_weight<-round(my_data$DSGN_WGT.y,0)
glm_model<-glm(EMPLOYED_208~.,data=my_data[,c("TREATMNT.y","FEMALE","AGE","BLACK","HISP","OTHERRAC",
                                                       "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                                                       "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                                                       "PERS_INC2","PERS_INC3","PERS_INC4",
                                                       
                                                       "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                                                       "WKEARNR","EMPLOYED_208"
)], family = binomial)



coefs<-round(coefficients(glm_model),6)
se<-round(coefficients(glm_model),6)
table2=round(summary(glm_model)$coefficients[,1:2],3)
write.csv(table2,paste0(my_path,"/JobCorps/STEP1_Replicate_Lee_2009/csv/table2.csv"))
