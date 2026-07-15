# Dietary attenuation
# 10–50% partial attenuation (the association is partly explained but a meaningful independent effect remains)
# >50% substantial attenuation (the adjustment variable explains a large share of the association)

bt <- function(x) paste0("`", x, "`")

results <- list()
failed <- list()

for (i in seq_len(nrow(candidate_paths))) {
  diet_var <- candidate_paths$diet[i]
  microbe_var <- candidate_paths$microbe[i]
  met_var <- candidate_paths$CHEM_ID[i]
  
  f1 <- as.formula(paste(bt(met_var), "~", bt(microbe_var),"+ covariates"))
  f2 <- as.formula(paste(bt(met_var), "~", bt(microbe_var), "+", bt(diet_var),"+ covariates"))
  m1 <- try(lm(f1, data = clr_s_met), silent = TRUE)
  m2 <- try(lm(f2, data = clr_s_met), silent = TRUE)
  
  if (inherits(m1, "try-error") | inherits(m2, "try-error")) {
    failed[[length(failed) + 1]] <- data.frame(
      row = i,
      diet = diet_label,
      diet_var = diet_var,
      microbe = microbe_var,
      metabolite = met_var,
      error_m1 = ifelse(inherits(m1, "try-error"), as.character(m1), NA),
      error_m2 = ifelse(inherits(m2, "try-error"), as.character(m2), NA))
    next}
  coef_m1 <- summary(m1)$coefficients
  coef_m2 <- summary(m2)$coefficients
  
  beta1 <- coef_m1[microbe_var, "Estimate"]
  beta2 <- coef_m2[microbe_var, "Estimate"]
  p1 <- coef_m1[microbe_var, "Pr(>|t|)"]
  p2 <- coef_m2[microbe_var, "Pr(>|t|)"]
  
  results[[length(results) + 1]] <- data.frame(
    diet = diet_label,
    diet_var = diet_var,
    microbe = microbe_var,
    metabolite = met_var,
    beta_m1 = beta1,
    beta_m2 = beta2,
    p_m1 = p1,
    p_m2 = p2,
    sig_m1 = p1 < 0.05,
    sig_m2 = p2 < 0.05,
    remains_significant = p1 < 0.05 & p2 < 0.05,
    attenuation = ((beta1 - beta2) / beta1) * 100,
    attenuation_abs = ((abs(beta1) - abs(beta2)) / abs(beta1)) * 100)}

attenuation_results <- bind_rows(results)
