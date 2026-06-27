# Predict secondary disease risk (Renal tubulo diseases) in diabetic patients
# using mediating proteins identified from the mediation analysis.
# Model: LASSO-Cox regression, evaluated by 10-fold cross-validation (C-index).

library(glmnet)
library(caret)
library(survival)

# --- Step 1: Select top mediating proteins ---
# Load mediation results and keep proteins with mediation ratio >= 10%
mediator_protein = read.csv("Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Result/Train test split significant protein.csv")
mediator_protein = subset(mediator_protein,ratio >= 0.1)


load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Data/protein clinical training id.RData")



# Focus on Diabetes -> Renal tubulo diseases pathway; rank by mediation ratio
selected_protein = subset(mediator_protein,X == "Diabetes" & Y == "RTID")
selected_protein = selected_protein[order(selected_protein$ratio,decreasing = T),]
Top10_protein = selected_protein$protein[1:5]

# --- Step 2: Extract protein expression for selected proteins ---
###Prediction by selected protein
protein_expression_selected = protein[,c("ID",Top10_protein)]
nmissing = apply(protein_expression_selected,2,function(x) return(sum(is.na(x))))
protein_expression_selected = na.omit(protein_expression_selected)




# --- Step 3: Load survival outcome data ---
# full_data contains survival time and event status for the Diabetes -> Renal tubulo diseases analysis
load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Data/Model 1/RTID as Y diabetes as x wide interval 9 confounders cleaned data.RData")
#load("~/Desktop/Postdoc/UKBioBank/UK biobank CVD MD ND mediation analysis/Wrap Up/Data/Model 1/HF as Y diabetes as x wide interval 9 confounders cleaned data.RData")



# Keep only diabetic patients (primary disease occurred) who also have protein data
full_data = subset(full_data,Occured == 1)
full_data = subset(full_data,ID %in% protein_expression_selected$ID)
sum(full_data$Occured_y == 2)  # number of renal disease events

# --- Step 4: Merge protein expression into survival data ---
full_data$ID = as.double(full_data$ID)
full_data = left_join(full_data,protein_expression_selected,by = "ID")
# Retain ID, survival time, event status, key covariates, and protein columns
full_data = full_data[,c(1,3,4,15:17,20:ncol(full_data))]
full_data = na.omit(full_data)

# --- Step 5: Prepare matrices for Cox model ---
# Y: survival outcome (time, status); X: covariates + protein expression
full_data$Occured_y = ifelse(full_data$Occured_y == 2,1,0)
full_data_training = subset(full_data,ID %in% training_id)
full_data_test = subset(full_data, !(ID %in% training_id))
Y_training = full_data_training[,c(3,2)]
Y_training = as.matrix(Y_training)
colnames(Y_training) = c("time","status")
X_training = full_data_training[,4:ncol(full_data_training)]
X_training = as.matrix(X_training)

Y_test = full_data_test[,c(3,2)]
Y_test = as.matrix(Y_test)
colnames(Y_test) = c("time","status")
X_test = full_data_test[,4:ncol(full_data_test)]
X_test = as.matrix(X_test)
model_0 = cv.glmnet(X_training,Y_training,family = "cox")
model_final = glmnet(X_training,Y_training,lambda = model_0$lambda.min,family = "cox")
risk = predict(model_final,X_test)
concordance <- survConcordance(Surv(Y_test) ~ risk)
C1_selected = concordance$concordance
C1_selected



