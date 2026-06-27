library(readr)
library(ggrepel)
library(ggplot2)
library(tidyverse)
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/Diabetes as X cleaned data wide interval.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/data/hypentensive 1 as X cleaned data wide interval.RData")

load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Data/protein clinical training id.RData")


full_data = as.data.frame(mood_disorder_x)
full_data = subset(full_data,Occured != 100)
table(full_data$Occured)
clinical$ID = as.character(clinical$ID)
full_data = left_join(full_data,clinical[,c(1,12,15:24)],by = "ID")


protein_list = colnames(protein)[-1]
#protein_list[2923]
p_value = numeric(0)
alpha = numeric(0)
protein$ID = as.character(protein$ID)
full_data = left_join(full_data,protein)


count = 1
gene = protein_list[1]
full_data = subset(full_data, ID %in% training_id)
X = full_data[,c(1,2,5:15)]

for (gene in protein_list) {
  y = as.data.frame(full_data[,gene])
  colnames(y) = "genei"
  model_data = cbind(X,y)
  model_data = na.omit(model_data)
  lmod = lm(genei ~ Occured + age_ins0+BMI_ins0+sex+white_ethnic+Sys+Dia+Smoking+Drinking+center_ins0+LDL,data = model_data)
  #lmod = lm(genei ~ Occured+Glucose + age_ins0+BMI_ins0+sex+Smoking+Drinking+LDL+center_ins0,data = model_data)
  a = summary(lmod)$coefficients
  p_value[count] = a[2,4]
  alpha[count] = a[2,1]
  count = count+1
}
sum(p_value <0.05)
p_adj = p.adjust(p_value,method = "bonferroni")
sum(p_adj<0.05)
rm(list=setdiff(ls(), c("p_value","alpha","p_adj","protein_list")))
