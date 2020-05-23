This archive contains scripts used in Semenova (2020) "Better Lee Bounds" to replicate Table "Estimate effects of PACES voucher  on test scores" (Table ABBKK)

The scripts are contained in sub‐directories Step1_Prepare_Data, Step2_Estimate_Bounds.

### Step 1: Basic data preparation
```prepare_data.R```  

-- takes ```aerdat4.sas7bdat``` (core dataset) and ```tab5v1.sas7bdat``` (test scores for test-takers only) from 
[Angrist Data Archive](https://economics.mit.edu/faculty/angrist/data1/data/angetal02)

-- merges them by ```ID```

-- adds indicator of missigness for each baseline covariate. Required R libraries are: ```tidyverse``` and ```sas7bdat'''

### Step 2: Estimate bounds 

```
module load R
Rscript estimate_bounds.R
```  

replicates Columns 2,3,4,5 from Table ABBKK

```Rscript estimate_lasso_bounds.R``` 

replicates Columns 6,7 from Table ABBKK. Output log is Table_ABBKK_2002.txt. Required R libraries are: ```hdm``` and ```expm```
