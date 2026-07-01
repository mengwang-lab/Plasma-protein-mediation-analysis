library(tidyverse)
library(survival)

# ============================================================
# File paths
# ============================================================
# Raw UK Biobank data are not included in this repository.
# The following files are placeholders for analysis-ready data.
# Users should generate these files after obtaining access to UK Biobank data.
# ============================================================

data_dir <- "data/processed"

x_file <- file.path(data_dir, "x_diabetes_cleaned_wide_interval.RData")

training_id_file <- file.path(data_dir, "protein_clinical_training_id.RData")

y_files <- list(
  cerebrovascular_disease = file.path(data_dir, "y_cerebrovascular_disease_cleaned_by_group.RData"),
  glomerular_diseases     = file.path(data_dir, "y_glomerular_diseases_cleaned_by_group.RData")
)

# Choose one outcome for analysis
selected_y <- "cerebrovascular_disease"
# selected_y <- "glomerular_diseases"


# ============================================================
# Load data
# ============================================================

load(x_file)
load(training_id_file)
load(y_files[[selected_y]])


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

