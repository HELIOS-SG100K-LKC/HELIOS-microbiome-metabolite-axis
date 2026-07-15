# Shannon diversity calculation
# Shannon assumes proportional data (each sample sums to 1) so normalise relative abundances per sample for shannon diversity calculation
library(vegan)
shannon_s <- data.frame(
  id = rel_s_fil$id,
  shannon_index = vegan::diversity(rel_p_fil2[,-c(1:3)]))

#---
# Partial correlation of Shannon index and metabolites
library(ppcor)

clr_s_met <- rel_s_fil2 %>%
  dplyr::select(-c(other)) %>%
  column_to_rownames("id") %>%
  clr() %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  inner_join(df) %>%
  inner_join(shannon_s)

covar <- clr_s_met %>%
  dplyr::select(covariates) %>%
  mutate(across(everything(), as.numeric))

partial_cor <- data.frame(metabolite = grep("^X", names(clr_s_met), value = TRUE),cor = NA_real_,p = NA_real_)

for (i in 1:nrow(partial_cor)) {
  result <- pcor.test(x = clr_s_met$shannon_index,
                      y = clr_s_met[[partial_cor$metabolite[i]]],
                      z = covar,
                      method = "spearman")
  partial_cor$cor[i] <- result$estimate
  partial_cor$p[i] <- result$p
}

partial_cor <- partial_cor %>%
  mutate(p_bon =p.adjust(p,method="bonferroni",n=length(p))) %>%
  filter(p_bon < 0.05)
