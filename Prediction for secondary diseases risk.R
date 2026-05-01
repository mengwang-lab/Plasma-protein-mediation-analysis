# Predict secondary disease risk (Renal tubulo diseases) in diabetic patients
# using mediating proteins identified from the mediation analysis.
# Model: LASSO-Cox regression, evaluated by 10-fold cross-validation (C-index).

library(glmnet)
library(caret)
library(survival)

# --- Step 1: Select top mediating proteins ---
# Load mediation results and keep proteins with mediation ratio >= 10%
hyper_dia = read.csv("Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Table/Hyper Diabetes Protein Mediation result.csv")
hyper_dia = subset(hyper_dia,ratio >= 0.1)

# Focus on Diabetes -> Renal tubulo diseases pathway; rank by mediation ratio
selected_protein = subset(hyper_dia,X == "Diabetes" & Y == "Renal tubulo diseases")
selected_protein = selected_protein[order(selected_protein$ratio,decreasing = T),]
Top10_protein = selected_protein$protein[1:10]

# --- Step 2: Extract protein expression for selected proteins ---
###Prediction by selected protein
load("~/Desktop/Postdoc/UKBioBank/wrap-up/Data/protein with confounders raw data.RData")
protein_expression_selected = protein_ins0[,c("ID",Top10_protein)]
nmissing = apply(protein_expression_selected,2,function(x) return(sum(is.na(x))))
protein_expression_selected = protein_expression_selected[,which(nmissing<10000)]
protein_expression_selected = na.omit(protein_expression_selected)




# --- Step 3: Load survival outcome data ---
# full_data contains survival time and event status for the Diabetes -> Renal tubulo diseases analysis
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Data/Renal tubulo diseases as Y diabetes as x wide interval 9 confounders cleaned data.RData")



# Keep only diabetic patients (primary disease occurred) who also have protein data
full_data = subset(full_data,Occured == 1)
full_data = subset(full_data,ID %in% protein_expression_selected$ID)
sum(full_data$Occured_y == 1)  # number of renal disease events

# --- Step 4: Merge protein expression into survival data ---
full_data$ID = as.double(full_data$ID)
full_data = left_join(full_data,protein_expression_selected,by = "ID")
# Retain ID, survival time, event status, key covariates, and protein columns
full_data = full_data[,c(1,3,4,9,12,13,19:ncol(full_data))]
full_data = na.omit(full_data)

# --- Step 5: Prepare matrices for Cox model ---
# Y: survival outcome (time, status); X: covariates + protein expression
Y = full_data[,c(3,2)]
Y = as.matrix(Y)
colnames(Y) = c("time","status")
X = full_data[,4:ncol(full_data)]
X = as.matrix(X)

# --- Step 6: 10-fold cross-validation with LASSO-Cox ---
# C-index measures discrimination ability (0.5 = random, 1 = perfect)
set.seed(234)
folds = createFolds(full_data$TT0_y,10)
C1_selected = numeric(0)
for (i in 1:10) {
  training = full_data[-folds[[i]],]
  test = full_data[folds[[i]],]
  training_x = training[,4:ncol(full_data)]
  training_x = as.matrix(training_x)
  training_y = training[,c(3,2)]
  training_y = as.matrix(training_y)
  colnames(training_y) = c("time","status")
  # Select optimal lambda by cross-validation, then fit final Cox-LASSO model
  model_0 = cv.glmnet(training_x,training_y,family = "cox")
  model_0$lambda.min
  modeli = glmnet(training_x,training_y,lambda = model_0$lambda.min,family = "cox")
  test_x = as.matrix(test[,4:ncol(full_data)])
  test_y = as.matrix(test[,3:2])
  colnames(test_y) = c("time","status")
  risk = predict(modeli,test_x)
  concordance <- survConcordance(Surv(test_y) ~ risk)
  C1_selected[i] = concordance$concordance
}
# Fold-level C-indices and overall mean
C1_selected
mean(C1_selected,na.rm = T)
