This archive contains scripts used in Semenova (2020) "Better Lee Bounds" to replicate Table "Estimate JobCorps effects on log hourly wage" (Table JobCorps)

The scripts are contained in sub‐directories Step1_Prepare_Data, Step2_Estimate_Bounds, Step3_Print_Figures

### Step 1: Basic data preparation
```prepare_main_data.R```  

-- takes ```key_vars.sas7bdat```, ```mileston.sas7bdat```, ```baseline.sas7bdat```, ```empl_tl.sas7bdat``` from 
[AER website](https://www.aeaweb.org/articles?id=10.1257/aer.98.5.1864)

-- merges them by ```MPRID```

-- recodes some covariates to match Lee(2009) and saves the result into data/dataLee2009.csv
Required R libraries are: ```tidyverse``` and ```sas7bdat```

```prepare_covariates.R```  
-- takes ```baseline.sas7bdat``` and ```dataLee2009.csv``` (see previous step)

-- identifies all categorical variables and recodes them as binary indicators. Output is not loaded at GitHub due to memory constraint.

### Step 2: Estimate bounds for Table JobCorps 

```
module load R
Rscript estimate_bounds.R
```  

replicates Columns 1,2,3,4 from Table JobCorps

```Rscript estimate_lasso_bounds.R``` 

replicates Columns 5,6 from Table JobCorps. Output log is Table_JobCorps.txt. 
Required R libraries are: ```hdm``` and ```expm```

### Step 3. Miscellaneous
```
Rscript evidence_against_monotonicity.R
```
computes estimates for Figure 1,2  
```
Rscript estimate_bounds_plb_nontreated_wage.R
```
computes estimates for Figure 3,5
### Step 4. Print Figures
```
Rscript Figure1.R Figure2.R Figure3.R Figure4.R Figure5.R 
```
