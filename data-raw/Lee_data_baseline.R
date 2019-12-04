my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data_baseline<-Lee_data[,c("MPRID","TREATMNT.y",baseline_varnames)]
devtools::use_data(Lee_data_baseline)
