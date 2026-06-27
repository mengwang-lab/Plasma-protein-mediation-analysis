# Install TwoSampleMR after activating the conda environment:
# conda activate manuscript-r-env
# Rscript install_TwoSampleMR.R

install.packages(
  "TwoSampleMR",
  repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org")
)

# Alternative if the command above fails, especially on some Linux/HPC systems:
# install.packages("remotes", repos = "https://cloud.r-project.org")
# remotes::install_github("MRCIEU/TwoSampleMR")
