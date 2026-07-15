# sCCA
library(PMA) 
library(ggtext)

clr_s_met <- rel_s_fil2 %>%
  dplyr::select(-c(other)) %>%
  column_to_rownames("id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  inner_join(df) 

X <- clr_s_met %>% dplyr::select(starts_with("X")) %>% as.matrix() # metabolites
X <- X[, colSums(!is.finite(X)) == 0, drop = FALSE] #
Y <- clr_s_met %>% dplyr::select(!starts_with("X")) %>% as.matrix() # microbes
dim(X)
dim(Y)
sum(!is.finite(X))
sum(!is.finite(Y))
X <- scale(X)
Y <- scale(Y)

# Using PMA::CCA, choose penaltyx/penaltyz between 0 and 1 (lower = more sparse)
set.seed(123)
perm_out <- PMA::CCA.permute(x = X,z = Y,typex = "standard",typez = "standard",nperms = 100)
best_penaltyx <- perm_out$bestpenaltyx
best_penaltyz <- perm_out$bestpenaltyz
scca_res <- PMA::CCA(x = X, z = Y, typex = "standard", typez = "standard", penaltyx = best_penaltyx, penaltyz = best_penaltyz, K = 1)  # K=1 canonical variate

# Non-zero coefficients (how strongly each selected feature contributes to the canonical axis; considered loadings)
# Top metabolites loadings
top_met <- tibble(feature = colnames(X),weight = as.numeric(scca_res$u)) %>%
  dplyr::filter(weight != 0) %>%
  dplyr::arrange(desc(abs(weight))) %>% 
  dplyr::mutate(direction = if_else(weight > 0, "Positive", "Negative"))
# Top microbes loadings
top_mic <- tibble(feature = colnames(Y),weight = as.numeric(scca_res$v)) %>%
  filter(weight != 0) %>%
  slice_max(abs(weight), n = top_n_mic) %>%
  arrange(desc(abs(weight))) %>%# 49 mic, 22 mic
  dplyr::mutate(direction = if_else(weight > 0, "Positive", "Negative")) %>%
  mutate(feature_short = sub("^([A-Za-z])[A-Za-z]*\\.(.*)$", "\\1.\\2", feature)) 

# Sample score calculation
sample_score <- tibble(id = rownames(X),met_score = as.numeric(X %*% scca_res$u),mic_score = as.numeric(Y %*% scca_res$v)) %>%
  mutate(joint_score = scale((mic_score + met_score) / 2)[,1],
         joint_cluster = if_else(joint_score >= 0, "Positive", "Negative"),
         joint_cluster = factor(joint_cluster, levels = c("Negative", "Positive")))
cor(sample_score$mic_score, sample_score$met_score, method = "pearson") # 0.7722886
cor_test <- cor.test(sample_score$mic_score, sample_score$met_score, method = "pearson") # 0.7722886
cor_test$estimate
cor_test$p.value
