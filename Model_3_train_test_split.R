# =============================================================================
# Mediation analysis: Diabetes (X) -> Protein (mediator) -> Cerebrovascular
# disease (Y), estimated on the TRAINING split of the UK Biobank cohort.
#
# For each candidate protein this script estimates:
#   - the DIRECT effect of diabetes on the outcome (Cox model), and
#   - the INDIRECT (mediated) effect through the protein (a * b product),
#     with a bootstrap distribution for inference.
#
# Designed to run as a SLURM array job: the full protein list is partitioned
# into chunks, one chunk per array task, and results are written per task.
# =============================================================================

# ---- Load inputs ------------------------------------------------------------
# Training-set sample IDs and the proteomic measurements.
load("~/UKBioBank/MA CVD MD/Training split mediation/Data/protein clinical training id.RData")
# Diabetes-associated proteins: provides `protein_list` and adjusted p-values `p_adj`.
load("~/UKBioBank/MA CVD MD/Training split mediation/Data/Diabetes associated proteins.RData")
# Outcome-specific analysis data (`full_data`). Switch the load below to change outcome:
#load("~/UKBioBank/MA CVD MD/Training split mediation/Data/GD as Y diabetes as x wide interval 9 confounders cleaned data.RData")
#load("~/UKBioBank/MA CVD MD/Training split mediation/Data/HF as Y diabetes as x wide interval 9 confounders cleaned data.RData")
load("~/UKBioBank/MA CVD MD/Training split mediation/Data/CereD as Y diabetes as x wide interval 9 confounders cleaned data.RData")

# ---- Array task index -------------------------------------------------------
#arr_id = 1
# Which chunk of proteins this job handles, read from the SLURM array index.
arr_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(arr_id)
#arr_id = 10
#arr_id = 1

library(survival)
library(tidyverse)

# ---- Select this chunk's proteins -------------------------------------------
# Keep only proteins significantly associated with diabetes (adjusted p < 0.05).
sub_protein = protein_list[which(p_adj<0.05)]
# Partition into 83-protein chunks; the last task (10) takes the remainder (748:827).
IDs= if (arr_id == 10) 748:827 else (83*(arr_id-1)+1):(83*(arr_id))
sub_protein = sub_protein[IDs]

# Attach this chunk's protein columns to the analysis data via sample ID.
protein$ID = as.character(protein$ID)
protein_sub = protein[,c("ID",sub_protein)]
full_data = left_join(full_data,protein_sub)

# ---- Bootstrap of the indirect (mediated) effect ----------------------------
# Resamples subjects 1000x and recomputes the a*b product each time, returning
# the bootstrap distribution of the indirect effect.
#   a = effect of diabetes (Occured) on the protein (linear model)
#   b = effect of the protein on the outcome, given diabetes (Cox model)
my_bootstrap_test = function(df){
  indirect_boot = numeric()
  for (k in 1:1000) {
    set.seed(k)                                  # fixed seed per replicate -> reproducible
    ids = sample(1:nrow(df),nrow(df),replace = T)# resample subjects with replacement
    dfi = df[ids,]
    # a-path: diabetes -> protein expression, adjusting for confounders.
    model2i = lm(genei ~ Occured + age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                   Dia+Sys+LDL,data = dfi)
    # b-path: protein -> outcome hazard, adjusting for diabetes and confounders.
    model3i = coxph(Surv(TT0_y,Occured_y)~Occured +genei+ age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                      Dia+Sys+LDL,data = dfi)
    ai = coefficients(model2i)[2]                # a = coefficient on Occured (diabetes)
    bi = coefficients(model3i)[2]                # b = coefficient on genei (protein)
    indirect_boot[k] = ai*bi                     # indirect effect = a * b
  }
  return(indirect_boot)
}

# ---- Accumulators for per-protein results -----------------------------------
count = 1
direct = numeric()     # point estimate of direct effect (Cox coef on diabetes)
direct_p = numeric()   # p-value of direct effect
indirect = numeric()   # bootstrap mean of the indirect effect
indirect_p = numeric() # bootstrap two-sided p-value of the indirect effect

# Restrict to training-set subjects, then keep only the columns used in models.
full_data = subset(full_data, ID %in% training_id)
X = full_data[,c(1,3:4,6,9:19)]

# ---- Loop over proteins -----------------------------------------------------
for (gene in sub_protein) {
  if (count%%10 == 1){
    print(count)                      # progress marker every 10 proteins
  }
  #print(gene)
  # Build a model frame: confounders/exposure (X) + this protein renamed to "genei".
  y = full_data[,gene]
  names(y) = "genei"
  model_data = cbind(X,y)
  colnames(model_data)[ncol(model_data)] = "genei"

  model_data = na.omit(model_data)    # complete-case analysis for this protein

  # Direct effect: diabetes -> outcome, controlling for the protein and confounders.
  cox = coxph(Surv(TT0_y,Occured_y)~Occured +genei+ age_ins0+BMI_ins0+white_ethnic+sex+Smoking+Drinking+
                Dia+Sys+LDL,data = model_data)
  result_cox = summary(cox)$coefficients
  direct[count] = result_cox[1,1]     # coefficient on Occured (diabetes)
  direct_p[count] = result_cox[1,5]   # its p-value

  # Indirect effect: bootstrap the a*b product.
  indirect_coefficient = my_bootstrap_test(model_data)
  indirect[count] = mean(indirect_coefficient)
  # Two-sided bootstrap p-value: proportion of replicates on the "wrong" side of 0.
  indirect_p[count] = 1- max(mean(indirect_coefficient>0),mean(indirect_coefficient<0))
  count = count+1
}

# ---- Collect and save results -----------------------------------------------
result = data.frame(protein = sub_protein,direct_effect = direct,direct_p = direct_p,
                    indirect = indirect,indirect_p = indirect_p)
rm(list=setdiff(ls(), c("arr_id","result")))   # drop everything except what we save
saveRDS(result,paste0("Diabetes as x CeD as y_training_set_20260515_",arr_id,".rds"))
