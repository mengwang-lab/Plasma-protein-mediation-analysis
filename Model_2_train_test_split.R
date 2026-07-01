library(readr)
library(ggrepel)
library(ggplot2)
library(tidyverse)

# ============================================================
# Load analysis-ready datasets
# ============================================================
# Raw UK Biobank data are not included in this repository.
# Users should prepare the required analysis-ready .RData files
# after obtaining access to UK Biobank data.
#
# Expected input files:
#   1. X exposure data, e.g., diabetes or hypertension
#   2. Training dataset containing participant IDs, clinical covariates,
#      and plasma protein expression data
# ============================================================


# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

data_dir <- "data/processed"

x_file <- file.path(data_dir, "x_diabetes_cleaned_wide_interval.RData")
# x_file <- file.path(data_dir, "x_hypertension_cleaned_wide_interval.RData")

protein_clinical_train_file <- file.path(
  data_dir,
  "protein_clinical_training_data.RData"
)


# ------------------------------------------------------------
# Load exposure data
# ------------------------------------------------------------

load(x_file)

# Example object expected after loading:
# x_data


# ------------------------------------------------------------
# Load training data with protein and clinical variables
# ------------------------------------------------------------

load(protein_clinical_train_file)

# Example object expected after loading:
# protein_clinical_train


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
