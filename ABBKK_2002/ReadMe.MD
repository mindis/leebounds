This archive contains scripts used in Semenova (2020) "Better Lee Bounds" to replicate Table "Estimate effects of PACES voucher 
on test scores".

The scripts are contained in sub‐directories Step1_Prepare_Data, Step2_Estimate_Bounds.

Step 1: Basic data preparation
'''prepare_data.R''' 
     --- takes the raw files from Angrist Data Archive 'https://economics.mit.edu/faculty/angrist/data1/data/angetal02'
         The files are '''aerdat4.sas7bdat''' (core dataset) and '''tab5v1.sas7bdat''' (test scores for test-takers only).
     --- merges the two files by 'ID'
     --- adds missingness indicators for each baseline covariate
     
Step 2: Estimate bounds
 
   
