# Metabolome variation
library(glmnet)  
library(caret)

clr_s_met <- rel_s_fil2 %>%
  dplyr::select(-c(other)) %>%
  column_to_rownames("id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  inner_join(df)

all_metabolites <- clr_s_met %>% 
  dplyr::select(X35:X999926119) # Check which metabolites have zero variance
constant_mets <- names(all_metabolites)[sapply(all_metabolites, function(x) var(x, na.rm = TRUE) == 0)]

microbe <- clr_s_met %>%
  dplyr::select(Klebsiella.pneumoniae:Faecalibacterium.sp..I2392)
x <- as.matrix(microbe)
covariate_matrix <- model.matrix(~ covariates - 1, data=clr_s_met)
x_full <-cbind(x,covariate_matrix)
y <- clr_s_met %>%
  dplyr::select(X35:X999926119) %>%
  dplyr::select(where(~ sum(!is.na(.)) > 1)) %>%
  dplyr::select(where(~ var(., na.rm = TRUE) > 0))

combined <- cbind(x_full, y)
combined_clean <- combined[complete.cases(combined), ]
x_clean <- as.matrix(combined_clean[, colnames(x_full)])
y_clean <- combined_clean[, colnames(y)]

results <- data.frame(metabolite = colnames(y_clean), lambda = NA, R2 = NA, mse = NA)

set.seed(123)
for (i in 1:ncol(y_clean)) {
  message(paste("Processing:", colnames(y_clean)[i]))
  y_response <- y_clean[[i]]
  
  tryCatch({
    ridge_model <- cv.glmnet(x_clean, y_response, alpha = 0,nfolds = 10, type.measure = "mse")
    lambda_1se <- ridge_model$lambda.1se
    model <- glmnet(x_clean, y_response, alpha = 0, lambda = lambda_1se)
    y_pred <- predict(model, x_clean,s = "lambda.1se")
    
    mse <- mean((y_response - y_pred)^2)
    R2_cv <- 1 - sum((y_response - y_pred)^2) / sum((y_response - mean(y_response))^2)
    
    results$lambda[i] <- lambda_1se
    results$R2_cv[i] <- R2_cv
    results$mse[i] <- mse
  }, error = function(e) {
    message(paste("Skipping metabolite:", colnames(y_clean)[i], " - Error:", e$message))
  })
}
results_merge <- results %>%
  rename(CHEM_ID_X=metabolite) %>%
  mutate(R2_cv_percent=R2_cv*100)

# Identify top metabolite per SUPER_PATHWAY
top_mets <- results_merge %>%
  group_by(SUPER_PATHWAY) %>%
  slice_max(order_by = R2_cv_percent, n = 1, with_ties = FALSE) %>%
  ungroup()
