# leebounds
This R package implements basic and sharp Lee (2009) bounds on treatments effects in the presence of selection/nonresponse bias. 


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
# References
David Lee "Training, wages, and sample selection: Estimating sharp bounds on treatment effects" Review of Economic Studies, 76(3):1071–1102.

Vira Semenova "Machine Learning for Set-Identified Linear Models" https://arxiv.org/abs/1712.10024 November 2019
# Support
Vira Semenova: semenovavira@gmail.com
