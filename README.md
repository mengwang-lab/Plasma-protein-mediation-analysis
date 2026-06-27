Cardiometabolic diseases such as diabetes, hypertension, and dyslipidemia often lead to serious complications affecting the heart, brain, kidneys, and lungs, yet the biological links between these conditions and their downstream effects are not well understood. This study examined whether proteins circulating in the blood help explain how these diseases progress and lead to further health problems. Using data from more than 50,000 individuals, the researchers identified hundreds of proteins that appear to connect primary diseases to later outcomes and found that some of these proteins may play a direct role in disease development. The findings suggest that blood-based proteins could help improve early detection and support more targeted approaches to preventing disease progression.

<img width="468" height="350" alt="image" src="https://github.com/user-attachments/assets/7db32ffa-3689-484a-9a41-19b339b81170" />


## Repository structure
````markdown
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
```

## License

MIT License.
