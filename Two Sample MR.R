# =============================================================================
# Configuration for two-sample MR analysis using UKB-PPP pQTL data and OpenGWAS
#
# Note:
# - Do not hard-code JWT tokens, passwords, or private paths in public scripts.
# - Set the OpenGWAS JWT as an environment variable instead.
# - Update the paths below according to your local or cluster environment.
# =============================================================================


# ---- OpenGWAS authentication ------------------------------------------------

# Recommended:
# Store your token in ~/.Renviron as:
# OPENGWAS_JWT=your_token_here
#
# Then restart R and load it using:
opengwas_jwt <- Sys.getenv("OPENGWAS_JWT")

if (opengwas_jwt == "") {
  stop("OpenGWAS JWT not found. Please set OPENGWAS_JWT in your environment.")
}


# ---- Project paths ----------------------------------------------------------

project_dir <- "path/to/project"

pqtl_path <- file.path(
  project_dir,
  "data/pqtl"
)

result_path <- file.path(
  project_dir,
  "results/two_sample_mr"
)

rsid_map_file <- file.path(
  project_dir,
  "data/reference/ref_olink_rsid_map_mac5_info03_b0_7_all_chr_patched_v2.csv"
)

plink_path <- file.path(
  project_dir,
  "tools/plink/plink"
)

bfile_path <- file.path(
  project_dir,
  "data/reference/ref_openGWAS_1kg_v3_LD_EUR/EUR"
)


# Load RSID map
rsid_map <- fread(rsid_map_file, select = c("ID", "rsid"))

# Set up Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#load("~/UKBioBank/gwas/heart failure proteins.RData")
#protein <- readRDS("~/UKBioBank/gwas/Glo protein.rds")
#protein <- readRDS("~/UKBioBank/gwas/Renal protein.rds")
protein <- c( "BCAN", "IL1RL1",  "NEFL", "PVR","APOE", "TNFSF10")
proteins = list()
proteins[[1]] = protein
names(proteins)[1] = "ALZ"

# Define the GWAS outcomes for neurological and mental diseases and cardiovascular diseases
outcome_ids <-  c("finn-b-G6_ALZHEIMER",	#Alzheimer disease
                  "finn-b-G6_ALZHEIMER_INCLAVO",	#Alzheimer disease, including avohilmo
                  "finn-b-G6_AD_WIDE",	#Alzheimer's disease, wide definition
                  "ebi-a-GCST90027158",	#Alzheimer's disease
                  "ieu-b-5067")	#Alzheimer's disease



outcomes <- available_outcomes(opengwas_jwt = opengwas_jwt)
desired_outcomes <- outcomes %>% 
  dplyr::filter(id %in% outcome_ids)
outcome_traits <- setNames(desired_outcomes$trait, desired_outcomes$id)
outcome_ncase <- setNames(desired_outcomes$ncase, desired_outcomes$id)
outcome_ncontrol <- setNames(desired_outcomes$ncontrol, desired_outcomes$id)

# Function definitions
get_gene_position <- function(gene_symbol) {
  tryCatch({
    gene_info <- getBM(
      attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters = 'hgnc_symbol',
      values = gene_symbol,
      mart = ensembl
    )
    gene_info <- gene_info %>% dplyr::filter(grepl("^[0-9]+$|^X$", chromosome_name))
    if (nrow(gene_info) > 0) {
      return(gene_info[1, ])
    } else {
      return(NULL)
    }
  }, error = function(e) {
    flog.error("Error getting gene positions for %s: %s", gene_symbol, e$message)
    return(NULL)
  })
}

# function to load and fiter exposure pQTL data
load_and_filter_pqtl <- function(protein, protein_set_name, pqtl_path, rsid_map) {
  protein_folder <- tolower(protein_set_name)
  pqtl_file <- file.path(pqtl_path, paste0("pqtl_", protein, ".csv"))
  
  if (!file.exists(pqtl_file)) {
    flog.error("File not found: %s", pqtl_file)
    return(NULL)
  }
  
  # Load and process pqtl file
  pqtl_raw <- fread(pqtl_file, select = c("ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "N", "LOG10P", "INFO")) %>%
    dplyr::mutate(CHROM = as.character(CHROM),
                  CHROM = ifelse(CHROM == "23", "X", CHROM))
  
  # Define special case gene mappings
  special_cases <- list(
    "NTPROBNP" = "NPPB",
    "C19ORF12" = "C19orf12",
    "BAP18" = "C17orf49",
    "SARG" = "C1orf116",
    "PALM2" = "PALM2AKAP2",
    "MYLPF" = "MYL11",
    "MENT" = "C1orf56",
    "LEG1" = "C6orf58",
    "GPR15L" = "C10orf99",
    "WARS" = "WARS1",
    "CERT" = "CERT1"
  )
  
  gene_name <- ifelse(protein %in% names(special_cases), special_cases[[protein]], protein)
  gene_position <- get_gene_position(gene_name)
  
  if (is.null(gene_position)) {
    flog.warn("No gene position found for protein: %s", protein)
    return(NULL)
  }
  
  cis_start <- max(gene_position$start_position - 1e6, 1)
  cis_end <- gene_position$end_position + 1e6
  chromosome <- gene_position$chromosome_name
  
  flog.info("Gene position for %s: Chrom%s from %d to %d", protein, chromosome, gene_position$start_position, gene_position$end_position)
  flog.info("%s gene cis locus: %d to %d", protein, cis_start, cis_end)
  
  instruments <- pqtl_raw %>%
    mutate(pval = 10^(-LOG10P)) %>%
    dplyr::filter(
      pval < 5e-8 & CHROM == chromosome & 
        GENPOS >= cis_start & GENPOS <= cis_end
    )
  
  if (nrow(instruments) == 0) {
    flog.warn("No instrumental variants found in the cis locus for %s.", protein)
    return(NULL)
  }
  
  instruments <- merge(instruments, rsid_map, by = "ID")
  missing_rsids <- sum(is.na(instruments$rsid))
  if (missing_rsids > 0) {
    flog.warn("%d variants without RSIDs map for %s.", missing_rsids, protein)
    instruments <- instruments[!is.na(instruments$rsid), ]
  }
  return(instruments)
}

## function to format filtered exposure pQTL data
format_mr_data <- function(instruments) {
  tryCatch({
    formatted_data <- format_data(
      dat = instruments,
      type = "exposure",
      snp_col = "rsid",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "pval",
      samplesize_col = "N",
      chr_col = "CHROM",
      pos_col = "GENPOS",
      info_col = "INFO",
      log_pval = FALSE
    )
    return(formatted_data)
  }, error = function(e) {
    flog.error("Error formatting MR data: %s", e$message)
    return(NULL)
  })
}

# function to clump the instrumental variants with plink1.9
perform_ld_clumping <- function(formatted_data, plink_path, bfile_path, clump_r2) {
  tryCatch({
    if (nrow(formatted_data) == 0) {
      flog.warn("No variants to clump.")
      return(data.frame())
    }
    
    clumped_results <- ld_clump(
      dplyr::tibble(
        rsid = formatted_data$SNP, 
        pval = formatted_data$pval.exposure,
        id = formatted_data$id.exposure
      ),
      plink_bin = NULL,
      bfile = NULL,
      clump_kb = 1000,  #sufficient to capture LD within the cis-region.
      clump_r2 = clump_r2,   #sufficient to capture LD within the cis-region.
      pop = "EUR",
      opengwas_jwt = opengwas_jwt
    )
    
    if (nrow(clumped_results) == 0) {
      flog.warn("Clumping resulted in no variants.")
      return(data.frame())
    }
    return(clumped_results)
  }, error = function(e) {
    flog.error("Error in clumping step: %s", e$message)
    return(data.frame())
  })
}

# function to filter instruments based on F-statistic
calculate_f_statistic <- function(final_instruments) {
  if (!"F_statistic" %in% names(final_instruments) && 
      "beta.exposure" %in% names(final_instruments) && 
      "se.exposure" %in% names(final_instruments)) {
    final_instruments <- final_instruments %>%
      mutate(F_statistic = (beta.exposure / se.exposure)^2)
  }
  filtered_instruments <- final_instruments %>%
    dplyr::filter(F_statistic >= 10)
  
  if (nrow(filtered_instruments) == 0) {
    flog.warn("All instruments have F-statistics less than 10. Consider reviewing your instrument selection.")
  }
  return(filtered_instruments)
}


# Function to process exposure instrument
read_format_exposure <- function(protein, protein_set_name, pqtl_path, rsid_map) {
  instruments <- load_and_filter_pqtl(protein, protein_set_name, pqtl_path, rsid_map)
  if (is.null(instruments) || nrow(instruments) == 0) {
    flog.warn("No instruments found for %s.", protein)
    return(NULL)
  }
  
  # Convert instruments to data.frame if it is a data.table
  if (inherits(instruments, "data.table")) {
    instruments <- as.data.frame(instruments)
  }
  
  formatted_data <- format_mr_data(instruments)
  if (is.null(formatted_data) || nrow(formatted_data) == 0) {
    flog.warn("Formatted data is empty for %s.", protein)
    return(NULL)
  }
  
  return(formatted_data)
}


# Function to process exposure instrument
perform_clumping_and_f_stat <- function(formatted_data, protein, clump_r2) {
  clumped_results <- perform_ld_clumping(formatted_data, plink_path, bfile_path, clump_r2)
  if (is.null(clumped_results) || nrow(clumped_results) == 0) {
    flog.warn("Clumped results are empty for %s with clump_r2 = %f.", protein, clump_r2)
    return(NULL)
  }
  
  # Merge formatted data with clumped results
  final_instruments <- merge(formatted_data, clumped_results, by.x = "SNP", by.y = "rsid")
  if (nrow(final_instruments) == 0) {
    flog.warn("No instruments after merging for %s.", protein)
    return(NULL)
  }
  
  final_instruments <- calculate_f_statistic(final_instruments)
  if (nrow(final_instruments) == 0) {
    flog.warn("No instruments with F-statistic >= 10 for %s.", protein)
    return(NULL)
  }
  
  return(final_instruments)
}

# Function to extract and harmonize outcome data
extract_and_harmonize_outcome_data <- function(final_instruments, outcome_id, opengwas_jwt, plink_path, bfile_path) {
  outcome_data <- extract_outcome_data(
    snps = final_instruments$SNP,
    outcomes = outcome_id,
    opengwas_jwt = opengwas_jwt,
    proxies = TRUE,
    rsq = 0.8,
    align_alleles = 1,
    palindromes = 1,
    maf_threshold = 0.3
  )
  
  if (is.null(outcome_data) || nrow(outcome_data) == 0) {
    flog.warn("No outcome data retrieved for outcome %s", outcome_id)
    return(NULL)
  }
  
  outcome_data <- as.data.frame(outcome_data)
  
  # Harmonize data
  harmonized_data <- harmonise_data(
    exposure_dat = final_instruments,
    outcome_dat = outcome_data,
    action = 2
  )
  
  if (nrow(harmonized_data) == 0) {
    flog.warn("No harmonized data available after matching exposure and outcome data.")
    return(NULL)
  }
  
  # Calculate r.outcome with correction for ascertainment bias
  ncase <- outcome_ncase[[outcome_id]]
  ncontrol <- outcome_ncontrol[[outcome_id]]
  
  harmonized_data$r.outcome <- get_r_from_lor(
    lor = harmonized_data$beta.outcome,
    af = harmonized_data$eaf.outcome,
    ncase = ncase,
    ncontrol = ncontrol,
    prevalence = ncase/(ncase + ncontrol),
    model = "logit",
    correction = TRUE
  )
  return(harmonized_data)
}

# Function to perform MR analysis
perform_mr_analysis <- function(harmonized_data, protein, result_path, trait_name) {
  
  # Create a list of datasets: original and Steiger-filtered
  datasets <- list(
    original = harmonized_data
    #steiger_filtered = steiger_filtering(harmonized_data)
  )
  
  # Initialize a list to store results for both datasets
  results <- list()
  
  # Loop over the datasets
  for (data_name in names(datasets)) {
    data <- datasets[[data_name]]
    
    # Remove duplicated SNPs
    data <- data[!duplicated(data$SNP), ]
    
    # Check the number of SNPs
    num_snps <- nrow(data)
    if (num_snps < 1) {
      flog.warn("No SNPs available for analysis in %s dataset for %s and %s.", data_name, protein, trait_name)
      results[[data_name]] <- NULL
      next
    }
    flog.info("Performing MR analysis on %s dataset (%d SNPs) for %s and %s.", data_name, num_snps, protein, trait_name)
    
    # Adjust method list based on the number of SNPs
    method_list <- if (num_snps == 1) {
      "mr_wald_ratio"
    } else if (num_snps < 3) {
      c("mr_ivw", "mr_ivw_fe", "mr_ivw_mre")
    } else {
      c("mr_ivw", "mr_ivw_fe", "mr_ivw_mre", "mr_ivw_radial",
        "mr_simple_median", "mr_weighted_median","mr_penalised_weighted_median",
        "mr_simple_mode", "mr_simple_mode_nome",
        "mr_egger_regression", 
        "mr_weighted_mode", "mr_weighted_mode_nome")
    }
    
    # Perform MR analysis
    mr_results <- tryCatch({
      mr(data, method_list = method_list)
    }, error = function(e) {
      flog.warn("MR analysis failed in %s dataset for %s and %s: %s", data_name, protein, trait_name, e$message)
      return(NULL)
    })
    
    # Sensitivity analyses (only if more than one SNP)
    sensitivity_results <- if (!is.null(mr_results) && num_snps > 1) {
      list(
        heterogeneity = tryCatch(mr_heterogeneity(data), error = function(e) NULL),
        pleiotropy = tryCatch(mr_pleiotropy_test(data), error = function(e) NULL),
        single_snp = tryCatch(mr_singlesnp(data), error = function(e) NULL),
        leave_one_out = if(num_snps > 2) tryCatch(mr_leaveoneout(data), error = function(e) NULL) else NULL
      )
    } else {
      list(heterogeneity = NULL, pleiotropy = NULL, single_snp = NULL, leave_one_out = NULL)
    }
    
    # Calculate F-statistics and mean F-statistic
    f_statistics <- (data$beta.exposure / data$se.exposure)^2
    mean_f_statistic <- mean(f_statistics, na.rm = TRUE)
    
    # Steiger directionality test
    steiger <- tryCatch(directionality_test(data), error = function(e) NULL)
    
    # Combine results for this dataset
    results[[data_name]] <- list(
      mr_results = mr_results,
      steiger_direction = steiger,
      mean_f_statistic = mean_f_statistic,
      sensitivity = sensitivity_results
    )
  }
  
  # Save results
  #saveRDS(results, file = file.path(result_path, paste0("MR_results_", protein, "_", trait_name, ".rds")))
  
  return(results)
}

# Function to process outcome data and perform MR analysis
process_outcome <- function(final_instruments, outcome_id, protein, clump_r2, trait_name, result_path, results) {
  flog.info("Analyzing %s as exposure to %s with clump_r2 = %f", protein, outcome_id, clump_r2)
  
  harmonized_data <- tryCatch({
    extract_and_harmonize_outcome_data(final_instruments, outcome_id, opengwas_jwt, plink_path, bfile_path)
  }, error = function(e) {
    flog.error("Error in harmonizing data for %s and %s with clump_r2 = %f: %s", protein, trait_name, clump_r2, e$message)
    return(NULL)
  })
  
  if (is.null(harmonized_data) || nrow(harmonized_data) == 0) {
    flog.warn("No harmonized data for %s and %s with clump_r2 = %f.", protein, trait_name, clump_r2)
    return(results)
  }
  
  tryCatch({
    # Perform MR analysis
    result <- perform_mr_analysis(harmonized_data, protein, result_path, trait_name)
    results[[paste0(protein, "_", outcome_id, "_", trait_name, "_LD_", clump_r2)]] <- result
  }, error = function(e) {
    flog.error("Error in MR analysis for %s and %s with clump_r2 = %f: %s", protein, trait_name, clump_r2, e$message)
  })
  
  return(results)
}

# Progress bar setup
total_iterations <- sum(sapply(proteins, length)) * length(outcome_ids)
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = total_iterations,
  width = 60
)

# Main analysis loop
results <- list()
for (protein_set_name in names(proteins)) {
  for (protein in proteins[[protein_set_name]]) {
    print(protein)
    flog.info("Processing protein: %s", protein)
    
    tryCatch({
      formatted_data <- read_format_exposure(protein, protein_set_name, pqtl_path, rsid_map)
      
      if (is.null(formatted_data)) {
        next
      }
      
      # Perform LD clumping and MR analysis with clump_r2 = 0.1 and 0.01
      for (clump_r2 in c(0.1, 0.01)) { 
        # Process instruments and calculate F-statistics
        final_instruments <- perform_clumping_and_f_stat(formatted_data, protein, clump_r2)
        
        if (is.null(final_instruments)) {
          next
        }
        
        # Perform MR analysis for each outcome with the current clump_r2
        for (outcome_id in outcome_ids) {
          trait_name <- outcome_traits[outcome_id]
          results <- process_outcome(final_instruments, outcome_id, protein, clump_r2, trait_name, result_path, results)
        }
      }
      
    }, error = function(e) {
      flog.error("Error processing protein %s: %s", protein, e$message)
    })
    
    pb$tick()  # Move progress bar tick here after all processing is done
  }
}

saveRDS(results, file = file.path(result_path, "MR_results_ALZ_diabetes_multi_level.rds"))
