This archive contains scripts used in Semenova (2020) "Better Lee Bounds" to replicate Tables "Medicaid effect on healthcare utilization" 
(Table Utilization) and "Medicaid effect on self-reported health" (Table Health). 

The scripts are contained in sub‐directories Step1_Prepare_Data, Step2_Std_Trim, Step3_ML, Step4_Print_Tables

### Step 1: Basic data preparation
```prepare_data.R```  

-- takes ```oregonhie_stateprograms_vars.dta```, ```oregonhie_ed_vars.dta``` from 
[OHIE](http://www.nber.org/oregon/4.data.html)

-- merges them by ```person_id```

-- identifies main outcomes of interest and recodes some covariates to match Finkelstein et al(2012) 
Required R libraries are: ```tidyverse``` and ```foreign```

### Step 2: Estimate lower bound on Medicaid effect (ITT and LATE frameworks) using standard trimming method

```
module load R
Rscript estimate_lower_bound_std_trim.R
Rscript bootstrap_lower_bound_std_trim.R
```  

computes estimated lower bound and standard errors for Medicaid effect on various outcomes of interest. Saves output to csv/ 


### Step 3. Step 2: Estimate lower bound on Medicaid effect (ITT and LATE frameworks) based on additional covariates selected by 
machine learning

```
Rscript estimate_lower_bound_ml.R
Rscript bootstrap_lower_bound_ml.R
```
computes estimated lower bound and standard errors for Medicaid effect on various outcomes of interest. Saves output to csv/ 

### Step 4. Print Tables
