# Data preprocessing and count matrix preparation
library (tidyverse)https://github.com/HELIOS-SG100K-LKC/HELIOS-microbiome-metabolite-axis/tree/main/1_Microbial_variation

# p = phylum counts
# g = genus counts
# s = species counts

p <- p %>%
  mutate(other=rowSums(dplyr::select(.,unclassified,others),na.rm=TRUE)) %>%
  dplyr::select(-unclassified,-others) #
g <- g %>%
  mutate(other=rowSums(dplyr::select(.,unclassified,others),na.rm=TRUE)) %>%
  dplyr::select(-unclassified,-others) 
s <- s %>%
  mutate(other=rowSums(dplyr::select(.,unclassified,others),na.rm=TRUE)) %>%
  dplyr::select(-unclassified,-others) 

# Function to convert sum-normalisation to relative abundance
process_data <- function(input_df){
  id <- input_df$SAM.id
  df <- input_df %>%
    dplyr::select(-SAM.id) %>%
    mutate_all(as.numeric)
  rowSum <- df %>%
    rowSums()
  df_rel <- df %>%
    mutate(across(everything(), ~./rowSum*100)) 
  df_rel <- cbind(id,df_rel)
  return(df_rel)
}

rel_p <- process_data(p) 
rel_g <- process_data(g) 
rel_s <- process_data(s) 

# Filter dataset to microbes present in at least 0.1% abundance and in 10% of cohort
threshold <- nrow(rel_p_fil)*0.1 #10% of cohort
list_sel=list()
count=1
for (d in list(rel_p_fil,rel_g_fil,rel_s_fil)){
  z=colSums(d>0.1)
  sel_col=row.names(as.data.frame(z[z>=threshold])) 
  list_sel[[count]]<-sel_col
  count=count+1
}
rel_p_fil2 <- rel_p_fil[,list_sel[[1]]] 
rel_g_fil2 <- rel_g_fil[,list_sel[[2]]] 
rel_s_fil2 <- rel_s_fil[,list_sel[[3]]] 

#----
# Bray-curtis dissimilarity and Aitchison distance calculation
library(vegan)
library(ecodist)

# Columns 1:3 are metadata columns
rowSums(rel_s_fil2[,-c(1:3)]) |> summary()
rel_s_fil2_join_norm <- rel_g_fil2 %>%
  rowwise() %>%
  mutate(across(-c(1:3), ~ .x / sum(c_across(-c(1:3))))) %>%
  ungroup()

# Version 1 to calculate aitchison distance 
microbe_clr_s <- rel_s_fil2 %>% dplyr::select(-c(1:3)) %>% as.matrix() %>% clr()
dsim_s <- dist(microbe_clr_s, method = "euclidean")

# Version 2 to calculate bray curtis
dsim_s <- vegdist(rel_s_fil2[,-c(1:3)],method="bray",diag=TRUE,upper=TRUE)

# PCoA
dsim_s_pcoa <- cmdscale(dsim_s,k=3,eig=TRUE)
dsim_s_pcoa_df <- as.data.frame(dsim_s_pcoa$points) %>%
  setNames(paste0("pcoa", 1:3)) %>%
  rownames_to_column("id")
pcoa_var <- dsim_s_pcoa$eig / sum(dsim_s_pcoa$eig[dsim_s_pcoa$eig > 0]) * 100

#----
# PERMANOVA variance partitioning
adonis_vars <- c(paste0("PC", 1:x), "other_covariates")
adonis_formula <- as.formula(paste("dsim_s ~", paste(adonis_vars, collapse = " + ")))
adonis_raw <- adonis2(adonis_formula,data = df,permutations = 999,by = "terms")
adonis_s <- as.data.frame(adonis_raw)
adonis_s$Variables <- rownames(adonis_s)
rownames(adonis_s_3) <- NULL

adonis_s <- adonis_s %>%
  mutate(R2_=R2*100,
         fdr_p = p.adjust(`Pr(>F)`, method = "BH", n = length(`Pr(>F)`)),
         label=case_when(
           fdr_p>0.05 ~"",
           fdr_p>0.01 ~ "*",
           fdr_p>0.001 ~ "**",
           !is.na(fdr_p) ~ "***",
           TRUE ~ NA_character_)) %>%
  group_by(Variables) %>%
  summarise(R2_ = sum(R2_, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(R2_)) %>%
  mutate(Category = case_when(
    Variables == "Metabolome PCs" ~ "Molecular",
    Variables == "Host genetic PCs" ~ "Genetic",
    Variables %in% c("Diet PCs", "Nutrient PCs", "Energy intake", "Probiotics","Current smoker", "Cigarettes/day", "Weight change") ~ "Lifestyle",
    Variables == "Medication count" ~ "Medication",
    Variables %in% c("Sex", "Age", "BMI", "Education level", "Employment status","Working hours", "Housing type", "Income") ~ "Sociodemographic",
    Variables %in% c("Breastfed", "Breastfed duration") ~ "Early life",
    Variables == "Sequencing batch" ~ "Technical factor",TRUE ~ "Other"))

# PERMDISP assessment of homogeneity of dispersion
betadisper(dsim_s, df$x)
