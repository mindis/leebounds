This archive contains scripts used in Semenova (2020) "Better Lee Bounds" to replicate Table "Estimate effects of PACES voucher  on test scores" (Table ABBKK)

The scripts are contained in sub‐directories Step1_Prepare_Data, Step2_Estimate_Bounds.

### Step 1: Basic data preparation
```prepare_data.R```  

-- takes ```aerdat4.sas7bdat``` (core dataset) and ```tab5v1.sas7bdat``` (test scores for test-takers only) from 
[Angrist Data Archive](https://economics.mit.edu/faculty/angrist/data1/data/angetal02)

-- merges them by ```ID```

-- adds indicator of missigness for each baseline covariate

### Step 2: Estimate bounds. Output log is Table_ABBKK_2002.txt

```estimate_bounds.R```  

replicates Columns 2,3,4 from Table ABBKK

```estimate_lasso_bounds.R``` 

replicates Columns 5,6 from Table ABBKK
