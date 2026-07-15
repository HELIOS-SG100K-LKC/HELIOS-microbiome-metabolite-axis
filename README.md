# HELIOS Microbiome–Metabolite Axis

This repository contains the analysis code used in the manuscript:
> A microbial ecological axis links gut metabolism with cardiometabolic health

---

## Repository structure
```
├── 1_Microbial_variation/
├── 2_Alpha_diversity/
├── 3_Metabolome_variation/
├── 4_Microbe_metabolite_association/
├── 5_sCCA/
├── 6_Dietary_attenuation/
├── 7_Cardiometabolic_association/
├── 8_Causal_inference/
└── README.md
```

Scripts are organised according to the order of analyses presented in the manuscript.

---

## Analysis pipeline

### 1. Microbial variation
To quantify the interindividual microbial variation and estimate the proportion of variance estimated by host, environmental and technical factors.
Script for
- Data preprocessing and count matrix preparation
- Bray-curtis dissimilarity and Aitchison distance calculation
- PCoA analysis
- PERMANOVA variance partitioning
- PERMDISP assessment of homogeneity of dispersion

### 2. Alpha diversity
Assess microbial diversity and association with host metabolites
Script for
- Shannon diversity calculation
- Partial Spearman correlations with plasma metabolites

### 3. Metabolome variation
Estimate the proportion of variance in individual plasma metabolites explained by the gut microbiome
Script for
-	Ridge regression with 10-fold cross-validation for lambda optimisation
-	Estimation of metabolite-specific R2

### 4. Microbe-metabolite association
Identify metabolites associated with microbial species and microbial functional pathways
Script for
-	Multivariable linear regression

### 5. Sparse canonical correlation analysis (sCCA)
Multi-omics analysis to identify coordinated microbe-metabolite axis
Script for
-	Sparse feature selection
-	Canonical loading estimation
-	Sample score calculation

### 6. Dietary attenuation
Assess whether dietary exposures attenuate microbe-metabolite associations
Script for
-	Attenuation of effect size

### 7. Cardiometabolic association
Epidemiological analysis between microbes and metabolites with cardiometabolic phenotypes
Script for
-	Multivariable linear and logistic regression

### 8. Causal inference
Prioritise candidate causal microbial metabolites
Script for
-	Two-sample Mendelian Randomisation
-	Mediation analysis

> Note: Metabolite GWAS summary statistics and genetic analysis scripts are available in companion repository ** HELIOS-Metabolite-genetics**.

## Software Analyses were performed in R (v4.5.1) using packages described in the Methods section of the manuscript.
