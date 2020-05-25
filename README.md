# leebounds
This R package implements sharp bounds on treatments effects in the presence of selection/nonresponse bias in randomized control trials. It includes basic Lee bounds  [Lee (2009)](https://academic.oup.com/restud/article-abstract/76/3/1071/1590707) and sharp bounds based on machine learning [Semenova(2019)](https://arxiv.org/abs/1712.10024). See vignette ```leebounds``` for description of the package (available through ```vignette("leebounds")``` once package is installed), and the package manual for documentation of the package functions.

# Installation
```
install.packages("devtools")
library(devtools) 
install_github("vsemenova/leebounds")
library(leebounds)
```
# Example
```
# compute basic Lee (2009) bounds for ATE in week 208
leedata=data.frame(treat=JobCorps_baseline$TREATMNT.y,selection=JobCorps_employment$week_208,outcome=JobCorps_wages$week_208)
leebounds(leedata)
```
# Support
Vira Semenova: semenovavira@gmail.com
