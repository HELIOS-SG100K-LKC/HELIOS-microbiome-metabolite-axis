# Multivariable linear regression
library(stats)

clr_s_met <- rel_s_fil2 %>%
  dplyr::select(-c(other)) %>%
  column_to_rownames("id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  inner_join(df) 

res_list <- list()

# Loop over response variables (metabolites)
for (response_col in 128:ncol(clr_s_met)) {
  response_var <- names(clr_s_met)[response_col]
  
  # Loop over predictor variables (microbes)
  for (predictor_col in 7:127) {
    predictor_var <- names(clr_s_met)[predictor_col]
    
    # Use tryCatch to safely handle any modeling errors
    tryCatch({
      model <- lm(as.formula(
        paste(response_var, "~", predictor_var,
              "+ covariates")),data = clr_s_met)
      
      model_summary <- summary(model)
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Response_Variable = response_var,
        Predictor_Variable = predictor_var,
        Beta_Estimate = model_summary$coefficients[2, "Estimate"],
        Std_Error = model_summary$coefficients[2, "Std. Error"],
        p = model_summary$coefficients[2, "Pr(>|t|)"]
      )
    }, error = function(e) {
      message(sprintf("Model failed for: %s ~ %s | Error: %s", response_var, predictor_var, e$message))
    })
  }
}
res_s_met <- do.call(rbind, res_list)

res_s_met <- res_s_met %>%
  mutate(bon_p = p.adjust(p, method = "bonferroni")) %>%
  mutate(bon_label = case_when(
    bon_p > 0.05 ~ "",
    bon_p <= 0.05 & bon_p > 0.01 ~ "*",
    bon_p <= 0.01 & bon_p > 0.001 ~ "**",
    bon_p <= 0.001 ~ "***",
    TRUE ~ NA_character_)) 
