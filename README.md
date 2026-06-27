In this study, we investigated the role of clinically accessible plasma proteomics as mediators in the development pathways from primary to secondary diseases under the causal inference framework using data from approximately 50,000 UK Biobank participants. Across three primary diseases (diabetes, hypertension, and dyslipidemia) and 18 subsequent conditions spanning four clinical domains, we identified 1,461 significant mediation pathways involving 395 unique plasma proteins.
<img width="468" height="350" alt="image" src="https://github.com/user-attachments/assets/7db32ffa-3689-484a-9a41-19b339b81170" />



One small issue: because this section contains nested code blocks, GitHub may occasionally format it incorrectly when copied all at once. A safer version is below, using four backticks around the full section:

````markdown
## Repository structure

```text
your-manuscript-repo/
├── README.md
├── env/
│   ├── manuscript_r_environment.yml
│   └── install_TwoSampleMR.R
├── R/
│   ├── Model_1_train_test_split.R
│   ├── Model_2_train_test_split.R
│   ├── Model_3_train_test_split.R
│   ├── mendelian_randomization.R
│   └── prediction_model.R
├── data/
│   └── README.md
└── results/
    └── README.md

````

## Environment setup

The R environment used for the analyses can be recreated using the provided conda environment file.

```bash
conda env create -f env/manuscript_r_environment.yml
conda activate manuscript-r-env
