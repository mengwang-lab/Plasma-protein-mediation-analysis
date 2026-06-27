library(tidyverse)
library(survival)


##Diabetes as x
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Diabetes as X cleaned data wide interval.RData")
#mood_disorder_x

##Load training ID
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Data/protein clinical training id.RData")




###Different Y
###CVD

#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Orther heart disease as Y cleaned data by group.RData")
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/cerebrovascular disease as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Artery disease as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Ischaemic heart disease as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank cvd MD ND mediation analysis/data/atherosclerosis as Y cleaned data by group.cvdata")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank cvd MD ND mediation analysis/data/heart failure as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Value disorders as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Myocardial Disorders as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Cardiac arrest and arrhythmias as Y cleaned data by group.RData")



##GD
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Renal tubulo diseases as Y cleaned data by group.RData")
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Glomerular diseases as Y cleaned data by group.RData")

##ND
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Parkison's disease as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Systemic atrophies primarily affecting the central nervous system as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Episodic and paroxysmal disorders of the central nervous system as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Diseases of myoneural junction and muscle as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Demyelinating diseases of the central nervous system as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Nerve nerve root and plexus disorders as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Other disorders of the nervous system as Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Alzheimer's disease as Y cleaned data by group.RData")

##RD
#load("~/Desktop/Postdoc/UKBioBank/UK biobank cvd MD ND mediation analysis/data/Chronic lower respiratory diseases Y cleaned data by group.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank cvd MD ND mediation analysis/data/Other_respiratory_diseases_principally_affecting_the_interstitium as Y cleaned data by group.RData")


cleaned_cvd = as.data.frame(cleaned_GD)
Diabetes_data = mood_disorder_x
Diabetes_data = as.data.frame(Diabetes_data)
full_data = left_join(cleaned_cvd,Diabetes_data,by = "ID")
full_data = subset(full_data,Occured != 100)
full_data_control = subset(full_data, Occured == 0)
full_data_treatment = subset(full_data, Occured == 1 & TT0<TT0_y)
full_data = rbind(full_data_control,full_data_treatment)



end_data = "09/01/2023"
end_data = as.character(end_data)
end_data = as.Date(end_data,format = "%m/%d/%Y")
full_data_event = subset(full_data,Occured_y ==1)
full_data_Null = subset(full_data,Occured_y == 0)
ins0 = full_data_Null$Ins0
ins0 = as.character(ins0)
ins0_null = as.Date(ins0,format = "%Y-%m-%d")
diff_null = end_data-ins0_null
diff_null = as.numeric(diff_null)
full_data_Null$TT0_y = diff_null
full_data = rbind(full_data_event,full_data_Null)
full_data$TT0_y = as.numeric(full_data$TT0_y)
table(full_data$Occured_y)
mean(full_data$Occured_y == 1)
clinical$ID = as.character(clinical$ID)
full_data = left_join(full_data,clinical[,c(1,12,15:24)],by = "ID")
full_data$Occured_y = as.factor(full_data$Occured_y)
full_data$Occured = as.factor(full_data$Occured)


full_data$age_ins0 = round(full_data$age_ins0)
full_data$Occured_y = as.numeric(full_data$Occured_y)
full_data_training = subset(full_data, ID %in% training_id)
cox = coxph(Surv(TT0_y,Occured_y) ~ Occured + age_ins0+white_ethnic+BMI_ins0+sex+Drinking+Smoking+Sys+Dia+LDL,data = full_data_training)
summary(cox)






rm(list=setdiff(ls(), c("full_data","training_id")))
save.image("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap UP/Data/Model 1/GD as Y diabetes as x wide interval 9 confounders cleaned data.RData")
