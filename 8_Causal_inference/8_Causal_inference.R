# Causal inference

# Two-sample Mendelian randomisation
library(TwoSampleMR)

# Using diastolic blood pressure as phenotype example
# Filter metabolite → phenotype df for these metabolites
# res_pheno_met dataframe is derived from 7_Cardiometabolic_association/7_Cardiometabolic_association.R
met_sig_pheno <- res_pheno_met %>%
  filter(Response_Variable == "Dbp", bon_p < 0.05)

# Read metabolite-gwas file - this file contains all SNPs associated with each metabolite 
all_loci_all<-vroom("~/../all_loci_all.txt")

# Read external summary statistics
GCST90310295<-fread("~/../SBP/GCST90310295.h.tsv")

# Filter GWAS for metabolites of interest
clump_input_GCST90310295 <- all_loci_all %>%
  filter(CHEM_ID %in% met_sig_pheno$Predictor_Variable) %>%
  filter(P_value<5e-08) %>% 
  filter(rsids %in% GCST90310295$rsid) %>% 
  dplyr::select(SNP = Variant_ID, P = P_value, CHEM_ID)
length(unique(clump_input_GCST90310295$CHEM_ID)) #N=48
out_dir <- "~/../clump_input_files"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

met_ids <- unique(clump_input_GCST90310295$CHEM_ID)
message("Writing instrument files for ", length(met_ids), " metabolites")

# Loop over each CHEM_ID and write a file for each metabolite
for (met in met_ids) {
  df_met <- clump_input_GCST90310295 %>% filter(CHEM_ID == met) %>% dplyr::select(-CHEM_ID)
  file_name <- paste0(met, ".txt")
  write.table(df_met,
              file = file.path(out_dir, file_name),
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
}

# Clump files using plink2 in bash in terminal

# Read clumped files
out_dir <- "~/../output"
clump_files <-list.files(out_dir, pattern = "\\.clumps$", full.names = TRUE)
clump_results <- lapply(clump_files, fread)
names(clump_results) <- gsub("_clumped.clumps", "",basename(clump_files))

clump_df <- bind_rows(
  lapply(names(clump_results), function(met) {
    clump_results[[met]] %>%
      mutate(CHEM_ID = met)})) 

# Filter met-gwas file
all_loci_clump_fil <- all_loci_all %>%
  inner_join(clump_df, join_by(Variant_ID==ID, CHEM_ID==CHEM_ID)) 

# Take list of rsIDs to keep and filters a GWAS summary statistics file so that only those SNPs are retained elsewhere in terminal
# Read filtered file
GCST90310295<-fread("~/../SBP/GCST90310295_byclump.txt") 

GCST90310295_5e08 <- GCST90310295 %>%
  filter(rsid %in% all_loci_clump_fil$rsids) # 44

# Prepare exposure data
exposure_dat <- format_data(all_loci_clump_fil_5e08,
                            type = "exposure",
                            snp_col = "rsids", 
                            beta_col = "Beta",
                            se_col = "SE",
                            eaf_col = "Effect_Allele_Freq",
                            effect_allele_col = "Effect_Allele",
                            other_allele_col = "Non_Effect_Allele",
                            pval_col = "P_value",
                            phenotype_col = "CHEM_ID") 

# Prepare outcome data
GCST90310295 <- as.data.frame(GCST90310295_5e08)
GCST90310295 <- GCST90310295 %>%
  rename(
    Variant_ID = rsid,
    Beta = beta,
    SE = standard_error,
    Effect_Allele = effect_allele,
    Non_Effect_Allele = other_allele,
    Effect_Allele_Freq = effect_allele_frequency,
    P_value = p_value) %>%
  mutate(Phenotype = "DBP") 
outcome_dat <- format_data(GCST90310295,
                           type = "outcome",
                           snp_col = "Variant_ID",
                           beta_col = "Beta",
                           se_col = "SE",
                           eaf_col = "Effect_Allele_Freq",
                           effect_allele_col = "Effect_Allele",
                           other_allele_col = "Non_Effect_Allele",
                           pval_col = "P_value",
                           phenotype_col = "Phenotype")
data <- harmonise_data(exposure_dat, outcome_dat, 
                       action = 1 )

# Run MR
mr_GCST90310295 <- mr(data,method_list = c("mr_ivw","mr_wald_ratio"))
length(unique(mr_GCST90310295$exposure))
mr_GCST90310295 <- mr_GCST90310295 %>%
  mutate(pval_fdr =p.adjust(pval, method = "fdr")) %>%
  left_join(all_loci_clump_fil %>%
              dplyr::select(CHEM_ID, CHEMICAL_NAME) %>% distinct(),by = c("exposure" = "CHEM_ID")) %>%
  mutate(
    b_lower = b - 1.96 * se,
    b_upper = b + 1.96 * se,
    OR = exp(b),                # odds ratio
    OR_lower = exp(b - 1.96*se), # 95% CI lower
    OR_upper = exp(b + 1.96*se), # 95% CI upper
    sig_fdr = factor(ifelse(pval_fdr < 0.05, "FDR < 0.05", "NS"),levels = c("NS", "FDR < 0.05"))) %>%
  filter(pval<0.05) %>%
  filter(pval_fdr<0.05)

# For metabolites with >=3 genetic instruments for sensitivity analysis (Egger/median/mode need >=3)
# Filter significant exposures from primary MR (already FDR-filtered)
sig_exposures <- mr_GCST90310295 %>%
  filter(pval_fdr < 0.05) %>%
  pull(exposure) %>% unique()
data_sig <- data %>% filter(exposure %in% sig_exposures)

multi <- data_sig %>% filter(mr_keep) %>% count(exposure) %>% filter(n >= 3) %>% pull(exposure)
data_multi <- data_sig %>% filter(exposure %in% multi)
mr_sensitivity <- mr(data_multi,
                     method_list = c("mr_ivw","mr_egger_regression",
                                     "mr_weighted_median","mr_weighted_mode")) %>%
  left_join(all_loci_clump_fil %>% dplyr::select(CHEM_ID, CHEMICAL_NAME) %>% distinct(),by = c("exposure" = "CHEM_ID")) %>%
  mutate(b_lower = b - 1.96*se, b_upper = b + 1.96*se)

#----
# Mediation analysis
library(mediation)
library(dplyr)
library(purrr)

covariates <- c("covariates")
cov_str <- paste(covariates, collapse = " + ")

# map phenotype labels to how they appear as columns in clr_s_met, and their type
pheno_type <- c(BMI = "cont", HbA1C = "cont", Sbp = "cont", Dbp = "cont",
                Hypertension = "bin", T2D = "bin")

run_mediation <- function(microbe, metabolite, phenotype, data, cov_str, sims = 1000) {
  phenotype <- as.character(phenotype)
  ptype <- pheno_type[[phenotype]]
  if (is.null(ptype)) return(NULL)
  
  # backtick-quote names (metabolite IDs start with X; microbe names have dots)
  m  <- sprintf("`%s`", microbe)
  md <- sprintf("`%s`", metabolite)
  y  <- sprintf("`%s`", phenotype)
  
  # mediator model: metabolite ~ microbe + covariates (always linear)
  f_med <- as.formula(sprintf("%s ~ %s + %s", md, m, cov_str))
  # outcome model: phenotype ~ microbe + metabolite + covariates
  f_out <- as.formula(sprintf("%s ~ %s + %s + %s", y, m, md, cov_str))
  
  fit_med <- tryCatch(lm(f_med, data = data), error = function(e) NULL)
  fit_out <- tryCatch(
    if (ptype == "bin") glm(f_out, data = data, family = binomial)
    else lm(f_out, data = data),
    error = function(e) NULL)
  
  if (is.null(fit_med) || is.null(fit_out)) return(NULL)
  
  med <- tryCatch(
    mediate(fit_med, fit_out,
            treat = gsub("`", "", m), mediator = gsub("`", "", md),
            sims = sims),
    error = function(e) NULL)
  if (is.null(med)) return(NULL)
  
  tibble(
    microbe = microbe, metabolite = metabolite, phenotype = phenotype,
    ACME       = med$d0,        # average causal mediation effect (indirect)
    ACME_p     = med$d0.p,
    ADE        = med$z0,        # average direct effect
    total      = med$tau.coef,
    prop_med   = med$n0,        # proportion mediated
    prop_med_p = med$n0.p)
}

#Create a list of valid triplets: microbe, metabolite, phenotype
# (i.e., only keep combinations that exist across the data sources)
triplet_df <- microbe_met_pairs %>%
  rename(microbe = Predictor_Variable, metabolite = Response_Variable) %>%
  inner_join(res_pheno2_fil %>% rename(microbe = Predictor_Variable, phenotype = Response_Variable),
             by = "microbe") %>%
  inner_join(res_pheno_met %>% rename(metabolite = Predictor_Variable, phenotype = Response_Variable),
             by = c("metabolite", "phenotype"))

# de-duplicate triplets before running
triplets <- triplet_df %>%
  distinct(microbe, metabolite, phenotype) %>%
  mutate(phenotype = as.character(phenotype)) %>% 
  filter(phenotype %in% names(pheno_type))   # only phenotypes we can model

# run all triplets (this is the slow step; sims=1000 x N triplets)
med_results <- pmap_dfr(
  list(triplets$microbe, triplets$metabolite, triplets$phenotype),
  ~ run_mediation(..1, ..2, ..3, data = clr_s_met, cov_str = cov_str))

# FDR-correct the mediation p-values
med_results %>% count(phenotype)   # Hypertension and T2D should now have rows
med_results <- med_results %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH"),
         sig_mediation = ACME_fdr < 0.05) %>%
  filter(ACME_fdr < 0.05) 

# harmonise phenotype labels between mediation (Sbp/Dbp/HbA1C) and MR (SBP/DBP/HbA1c)
pheno_map <- c(BMI="BMI", HbA1C="HbA1c", Sbp="SBP", Dbp="DBP",
               Hypertension="Hypertension", T2D="T2D")

med_tiered <- med_results %>%
  mutate(phenotype = as.character(phenotype),
         phenotype_mr = pheno_map[phenotype]) 
