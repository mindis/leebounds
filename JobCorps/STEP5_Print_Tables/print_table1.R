rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/leebounds/"

source(paste0(my_path,"/R/utils.R"))

estimates<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_estimates.csv"),row.names=1)
sd<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_CR.csv"),row.names=1)
im<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_CR_IM.csv"),row.names=1)

table4latex<-print_table(estimates,sd,im,digs=3)
write.table(print(xtable(table4latex, type="latex"),include.rownames =FALSE ),paste0(my_path,"/JobCorps/STEP5_Print_Tables/Table1.txt"),append=FALSE)

estimates<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_estimates.csv"),row.names=1)
sd<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_CR.csv"),row.names=1)
im<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_IM.csv"),row.names=1)

estimates_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_estimates_lasso.csv"),row.names=1)
sd_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_CR_lasso.csv"),row.names=1)
im_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col1235_IM_lasso.csv"),row.names=1)


estimates<-cbind(estimates[,1:6], estimates_lasso[,1:2],estimates[,7:8])
sd<-cbind(sd[,1:6], sd_lasso[,1:2],sd[,7:8])
im<-cbind(im[,1:6], im_lasso[,1:2],im[,7:8])



table4latex<-print_table(estimates,sd,im,digs=3)

write.table(print(xtable(table4latex), type="latex",include.rownames =FALSE ),paste0(my_path,"/JobCorps/STEP5_Print_Tables/Table2.txt"),append=FALSE)