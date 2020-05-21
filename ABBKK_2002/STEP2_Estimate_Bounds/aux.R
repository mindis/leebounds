report_results<-function(estimates,CR) {
  estimates<-lapply(estimates,round,3)
  CR<-lapply(CR,round,3)
  ### Panel A: MATH
  print(paste0("Panel A: MATH, Estimate: [", estimates$MATH[1]," ", estimates$MATH[2],"]" ))
  print(paste0("Panel A: MATH, CR: (", CR$MATH[1]," ", CR$MATH[2],")" ))
  ### Panel B: READING
  print(paste0("Panel B: READING, Estimate: [", estimates$READING[1]," ", estimates$READING[2],"]" ))
  print(paste0("Panel B: READING  CR: (", CR$READING[1]," ", CR$READING[2],")" ))
  
  ### Panel C: WRITING
  print(paste0("Panel C: WRITING,  Estimate: [", estimates$WRITING[1]," ", estimates$WRITING[2],"]" ))
  print(paste0("Panel C: WRITING   CR: (", CR$WRITING[1]," ", CR$WRITING[2],")" ))
  
}