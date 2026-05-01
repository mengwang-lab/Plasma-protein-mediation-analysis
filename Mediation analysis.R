# ============================================================
# Mediation Analysis (example): Diabetes (X) → Protein → Heart Failure (Y)
# Model 3: Tests whether protein expression mediates the effect
# of diabetes on CVD outcomes, using SLURM array parallelization.
# Each array job processes a subset of significant proteins.
# ============================================================

# Load pre-filtered data from Model 2 (only the active outcome is uncommented)
load("~/UKBioBank/MA CVD MD/Data/9 confounders/model 2 result diabetes as x heart failure as y.RData")



# Read the SLURM array task ID to determine which protein chunk this job handles
arr_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(arr_id)
#arr_id = 10
#arr_id = 1
library(survival)
# Restrict to proteins that passed multiple-testing correction in Model 2
sub_protein = protein_list[which(p_adj<0.05)]
# Convert categorical covariates to factors for correct model encoding
full_data$center_ins0 = as.factor(full_data$center_ins0)
full_data$Smoking = as.factor(full_data$Smoking)
full_data$Drinking = as.factor(full_data$Drinking)

# -------------------------------------------------------
# Bootstrap function for the indirect (mediation) effect
# Estimates the product-of-coefficients (a*b) via 1000
# bootstrap resamples of the input data frame.
#   a = effect of Diabetes (Occured) on protein expression (linear model)
#   b = effect of protein on heart failure risk adjusting for Diabetes (Cox)
#   indirect = a * b (the mediated portion of the total effect)
# -------------------------------------------------------
my_bootstrap_test = function(df){
  indirect_boot = numeric()
  for (k in 1:1000) {
    set.seed(k)
    # Resample with replacement to get bootstrap dataset
    ids = sample(1:nrow(df),nrow(df),replace = T)
    dfi = df[ids,]
    # Path a: linear model — how much does Diabetes change the protein level?
    model2i = lm(genei ~ Occured + age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                   Dia+Sys,data = dfi)
    # Path b: Cox model — does the protein predict heart failure after adjusting for Diabetes?
    model3i = coxph(Surv(TT0_y,Occured_y)~Occured +genei+ age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                      Dia+Sys,data = dfi)
    ai = coefficients(model2i)[2]   # coefficient for Occured in path-a model
    bi = coefficients(model3i)[2]   # coefficient for genei in path-b model
    indirect_boot[k] = ai*bi        # indirect effect estimate for this bootstrap sample
  }
  return(indirect_boot)
}

# -------------------------------------------------------
# Initialise result containers and define the protein chunk
# for this SLURM array job (126 proteins per job, last job
# handles the remainder up to index 1260).
# -------------------------------------------------------
count = 1
direct = numeric()
direct_p = numeric()
indirect = numeric()
indirect_p = numeric()
# Covariate columns: time-to-event, event flag, age, BMI, ethnicity, sex,
# smoking, drinking, diastolic BP, systolic BP, and two disease indicators
X = full_data[,c(3:4,6,9:18)]
IDs= if (arr_id == 10) 1135:1260 else (126*(arr_id-1)+1):(126*(arr_id))
sub_protein = sub_protein[IDs]

# -------------------------------------------------------
# Main loop: for each protein, estimate the direct effect
# (Diabetes → Heart failure, protein included as covariate)
# and the indirect effect (via bootstrapped mediation).
# -------------------------------------------------------
for (gene in sub_protein) {
  if (count%%10 == 1){
    print(count)   # progress tracker: prints every 10 proteins
  }
  #print(gene)
  y = full_data[,gene]
  names(y) = "genei"
  model_data = cbind(X,y)
  colnames(model_data)[ncol(model_data)] = "genei"
  model_data = na.omit(model_data)   # drop rows with missing values

  # Full Cox model: direct effect of Diabetes on heart failure with protein as covariate
  cox = coxph(Surv(TT0_y,Occured_y)~Occured +genei+ age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                Dia+Sys,data = model_data)
  result_cox = summary(cox)$coefficients
  direct[count] = result_cox[1,1]    # log hazard ratio for Diabetes (direct effect)
  direct_p[count] = result_cox[1,5]  # p-value for the direct effect

  # Bootstrap-based indirect effect and its empirical p-value
  indirect_coefficient = my_bootstrap_test(model_data)
  indirect[count] = mean(indirect_coefficient)
  # Two-sided p-value: 2 * min(P(indirect>0), P(indirect<0))
  indirect_p[count] = 1- max(mean(indirect_coefficient>0),mean(indirect_coefficient<0))
  count = count+1
}

# Collect results and save — clean up workspace to free memory before saving
result = data.frame(protein = sub_protein,direct_effect = direct,direct_p = direct_p,
                    indirect = indirect,indirect_p = indirect_p)
rm(list=setdiff(ls(), c("arr_id","result")))
saveRDS(result,paste0("Diabetes as x heart failure as y_20241128_",arr_id,".rds"))
