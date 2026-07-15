# Cardiometabolic association
library(lme4)
library(car)

# Sig microbe-metabolite pairs
sig_pairs <- res_s_met %>% 
  filter(bon_p<0.05)
sig_met <-unique(sig_pairs$CHEMICAL_NAME) 
sig_microbe <-unique(sig_pairs$Predictor_Variable) 

# Filter significant microbes to assess cardiometabolic association
clr_s_met <- rel_s_fil2 %>%
  dplyr::select(-c(other)) %>%
  column_to_rownames("id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(df)

clr_s_met_fil<-clr_s_met %>%
  dplyr::select(id,covariates,phenotypes,
                all_of(unique(sig_pairs$Predictor_Variable))) %>%
  mutate(across(c(Hypertension,T2D), as.factor)) %>%
  mutate_at(c("HbA1C","BMI","Sbp","Dbp"),~(scale(.) %>% as.vector)) 

col_names <- c("Response_Variable", "Predictor_Variable", "Beta_Estimate", "Std_Error", "p")
res_pheno <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(res_pheno) <- col_names

# Linear & logistic regression models
total <- length(names(clr_s_met_fil)[6:11]) * length(unique(sig_pairs$Predictor_Variable))
pb <- txtProgressBar(min = 0, max = total, style = 3)
counter <- 0

for (response_var in names(clr_s_met_fil)[6:11]) {
  y <- clr_s_met_fil[[response_var]]
  is_binary <- is.factor(y) && length(levels(y)) == 2 || all(y %in% c(0, 1), na.rm = TRUE)
  
  for (predictor_var in unique(sig_pairs$Predictor_Variable)) {
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
    
    formula_str <- paste0(response_var, " ~ ", predictor_var, " + covariates")
    
    model <- tryCatch({
      if (is_binary) {
        glm(as.formula(formula_str), data = clr_s_met_fil, family = binomial)
      } else {
        lm(as.formula(formula_str), data = clr_s_met_fil)
      }
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      model_summary <- summary(model)
      coefs <- model_summary$coefficients
      
      # Check if predictor_var is present in coefficients
      if (predictor_var %in% rownames(coefs)) {
        beta_estimate <- coefs[predictor_var, "Estimate"]
        std_error <- coefs[predictor_var, "Std. Error"]
        p <- if (is_binary) coefs[predictor_var, "Pr(>|z|)"] else coefs[predictor_var, "Pr(>|t|)"]
        
        new_row <- data.frame(
          Response_Variable = response_var,
          Predictor_Variable = predictor_var,
          Beta_Estimate = beta_estimate,
          Std_Error = std_error,
          p = p
        )
        
        res_pheno <- rbind(res_pheno, new_row)
      }
    }
  }
}
res_pheno<-res_pheno %>%
  mutate(fdr_p = p.adjust(p, method = "BH"),
         bon_p = p.adjust(p, method = "bonferroni")) %>%
  filter(fdr_p<0.05) %>% 
  mutate(log_fdr_p = -log10(fdr_p),
         Predictor_Variable = factor(Predictor_Variable),
         Response_Variable = factor(Response_Variable))

# Filter significant metabolites to assess cardiometabolic association
clr_s_met_fil<-clr_s_met %>%
  dplyr::select(id,covariates,phenotypes,
                all_of(unique(sig_pairs$Response_Variable))) %>%
  filter(FREG5_Ethnic_Group!="O") %>%
  mutate_at(c("HbA1C","BMI","Sbp","Dbp"),~(scale(.) %>% as.vector)) 

col_names <- c("Response_Variable", "Predictor_Variable", "Beta_Estimate", "Std_Error", "p")
res_pheno_met <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(res_pheno_met) <- col_names

# Linear & logistic regression models
total <- length(names(clr_s_met)[5:10]) * length(unique(sig_pairs$Response_Variable))
pb <- txtProgressBar(min = 0, max = total, style = 3)
counter <- 0

for (response_var in names(clr_s_met)[5:10]) {#
  y <- clr_s_met[[response_var]]
  is_binary <- (is.factor(y) && length(levels(y)) == 2)
  
  for (predictor_var in unique(sig_pairs$Response_Variable)) {
    #if (predictor_var == response_var) next 
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
    
    formula_str <- paste0("`", response_var, "` ~ `", predictor_var, "` + FREG5_Ethnic_Group + FREG7_Gender + FREG8_Age + SSID")
    
    model <- tryCatch({
      if (is_binary) {
        glm(as.formula(formula_str), data = clr_s_met, family = binomial)
      } else {
        lm(as.formula(formula_str), data = clr_s_met)
      }
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      model_summary <- summary(model)
      coefs <- model_summary$coefficients
      
      rows <- grep(paste0("^`?", predictor_var), rownames(coefs))
      
      if (length(rows) > 0) {
        for (r in rows) {
          beta_estimate <- coefs[r, "Estimate"]
          std_error <- coefs[r, "Std. Error"]
          p <- if (is_binary) coefs[r, "Pr(>|z|)"] else coefs[r, "Pr(>|t|)"]
          
          new_row <- data.frame(
            Response_Variable = response_var,
            Predictor_Variable = rownames(coefs)[r],
            Beta_Estimate = beta_estimate,
            Std_Error = std_error,
            p = p
          )
          
          res_pheno_met <- rbind(res_pheno_met, new_row)
        }
      }
    }
  }
}
res_pheno_met<-res_pheno_met %>%
  mutate(fdr_p = p.adjust(p, method = "BH"),
         bon_p = p.adjust(p, method = "bonferroni")) %>%
  filter(bon_p<0.05) %>% 
  mutate(log_fdr_p = -log10(fdr_p),
         Predictor_Variable = factor(Predictor_Variable),
         Response_Variable = factor(Response_Variable))
