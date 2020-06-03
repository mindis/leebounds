rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")
library(reldist)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
my_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
my_data<-as.data.frame(my_data)

source(paste0(my_path,"/JobCorps/STEP1_Replicate_Lee_2009/utils.R"))

my_data_control<-filter(my_data,TREATMNT.y==0)
my_data_treat<-filter(my_data,TREATMNT.y==1)
sink(paste0(my_path,"/JobCorps/STEP1_Replicate_Lee_2009/Table4.log"))
print("Control: Number of obs")
print(dim(my_data_control)[1])
print("Control: Prop non missing")
print(weighted.mean(my_data_control$EARNH208>0,w=my_data_control$DSGN_WGT.y))
print("Control: Mean log wage")
log_wage_control<-log(my_data_control$EARNH208[my_data_control$EARNH208>0 & my_data_control$HWH208>0 ]/my_data_control$HWH208[my_data_control$EARNH208>0 & my_data_control$HWH208>0])
print(weighted.mean(log_wage_control,
                    w=my_data_control$DSGN_WGT.y[my_data_control$EARNH208>0 & my_data_control$HWH208>0]))

print("Treat: Number of obs")
print(dim(my_data_treat)[1])
print("Treat: Prop non missing")
print(weighted.mean(my_data_treat$EARNH208>0,w=my_data_treat$DSGN_WGT.y))
print("Treat: Mean log wage")
log_wage_treat<-log(my_data_treat$EARNH208[my_data_treat$EARNH208>0]/my_data_treat$HWH208[my_data_treat$EARNH208>0])
print(weighted.mean(log_wage_treat,
                    w=my_data_treat$DSGN_WGT.y[my_data_treat$EARNH208>0]))

p=weighted.mean(my_data_control$EARNH208>0,w=my_data_control$DSGN_WGT.y)/weighted.mean(my_data_treat$EARNH208>0,w=my_data_treat$DSGN_WGT.y)

print("Trimming threshold")
print(round(p,3))

### Table 6 ###
prop<-rep(0,4)
prop_treat<-rep(0,4)
prop_control<-rep(0,4)
bounds<-matrix(0,nrow=2,ncol=4)
untrimmed_effect<-rep(0,4)
selected_weeks<-c(45,90,135,180)

for (i in 1:length(selected_weeks)) {
  # prepare data
  week<-selected_weeks[i]
  print (paste0("Results for week ", week))
  
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  
  logwage_week<-as.numeric(as.matrix(log(my_data[,earn_name]/my_data[,hwh_name])))

  
  leedata_week<-as.data.frame(cbind(treat=my_data$TREATMNT.y,
                                    selection=as.numeric(my_data[,earn_name]>0),
                                    outcome=logwage_week,
                                    weights=my_data$DSGN_WGT.y))

  leedata_week_control<-filter( leedata_week,treat==0)
  leedata_week_treat<-filter( leedata_week,treat==1)
  
  prop_control[i]<-weighted.mean( leedata_week_control$selection==1,w= leedata_week_control$weights)
  prop_treat[i]<-weighted.mean( leedata_week_treat$selection==1,w= leedata_week_treat$weights)
  
  res<-leebounds_unknown_sign(leedata_week,treat_helps=(week>=90))
  
  
  prop[i]<-GetThresh(res)
  bounds[,i]<-GetBounds(res)
  untrimmed_effect[i]<-  weighted.mean( logwage_week[my_data$TREATMNT.y==1 & leedata_week$selection==1],w=my_data$DSGN_WGT.y[my_data$TREATMNT.y==1 & leedata_week$selection==1])-
    weighted.mean( logwage_week[my_data$TREATMNT.y==0 & leedata_week$selection==1],w=my_data$DSGN_WGT.y[my_data$TREATMNT.y==0 & leedata_week$selection==1])
}
bounds<-apply(t(bounds),2,round,3)
colnames(bounds)<-c("lower","upper")
rownames(bounds)<-selected_weeks
selected_weeks
sink(file=NULL)