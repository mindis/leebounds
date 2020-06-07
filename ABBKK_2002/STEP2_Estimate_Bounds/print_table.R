rm(list=ls())
## load relevant libraries
library(xtable)
## specify path to code 
my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/utils.R"))

load(paste0(my_path,"/ABBKK_2002/draft_older/ABBKK_2002_baseline.RData"))
load(paste0(my_path,"/ABBKK_2002/draft_older/ABBKK_2002_Col67.RData"))

angrist_estimates<-matrix(c(0.178,0.178, 0.204,0.204,0.126, 0.126), ncol=2,byrow=TRUE)
angrist_cr<-matrix(c(-0.058, 0.413, -0.021, 0.429, -0.101, 0.353), ncol=2,byrow=TRUE)
estimates_table<-cbind(angrist_estimates,
                       matrix(unlist(estimates),ncol=2,byrow=TRUE),
                       matrix(unlist(estimates_nonmonotone),ncol=2,byrow=TRUE),
                       matrix(unlist(orthoestimates),ncol=2,byrow=TRUE),
                       matrix(unlist(orthoestimates2),ncol=2,byrow=TRUE),
                       matrix(unlist(orthoestimates_postlasso),ncol=2,byrow=TRUE),
                       matrix(unlist(estimates_plb),ncol=2,byrow=TRUE)
                       )


colnames(estimates_table)<-c("ABBKK_lb","ABBKK_ub",
                             "Lee_2009_lb","Lee_2009_ub",
                             "Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone",
                             "OrthoLee_28_lb","OrthoLee_28_ub",
                             "OrthoLee_28_lb2","OrthoLee_28_ub2",
                             "OrthoLee_28_postlasso","OrthoLee_28_postlasso",
                             "PLB_LB","PLB_UB")

CR_table<-cbind(angrist_cr,
  matrix(unlist(CR),ncol=2,byrow=TRUE),
  matrix(unlist(CR_nonmonotone),ncol=2,byrow=TRUE),
                       matrix(unlist(CR_ortho),ncol=2,byrow=TRUE),
                       matrix(unlist(CR_ortho2),ncol=2,byrow=TRUE),
                       matrix(unlist(CR_ortho_postlasso),ncol=2,byrow=TRUE),
                       matrix(unlist(CR_plb),ncol=2,byrow=TRUE)
)
colnames(CR_table)<-c("ABBKK_lb","ABBKK_ub",
                      "Lee_2009_lb","Lee_2009_ub",
                      "Lee_2009_lb_nonmonotone","Lee_2009_ub_nonmonotone",
                      "OrthoLee_28_lb","OrthoLee_28_ub",
                      "OrthoLee_28_lb2","OrthoLee_28_ub2",
                      "OrthoLee_28_postlasso","OrthoLee_28_postlasso",
                      "PLB_LB","PLB_UB")




estimates_table<-apply(estimates_table,2,round,3)
CR_table<-apply(CR_table,2,round,3)
table4latex<-print_table(estimates_table,CR_table,digs=3)

write.table(print(xtable(table4latex), type="latex",include.rownames =FALSE ),paste0(my_path,"/ABBKK_2002/STEP2_Estimate_Bounds/Table1.txt"),append=FALSE)