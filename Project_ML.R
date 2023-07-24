#library(DESeq2)
library(ggplot2)
library(tidyverse)
library(tidymodels)
library(Seurat)
#library(scCustomize)
#library(airway)
library(limma)
library(pheatmap)
library(ggsci)
library(genefilter)
library(dplyr)
library(randomForest)
library(stats)
library(tibble)
library(ggfortify)
library(affy)
library(ROCR)
library(Hmisc)
library(patchwork)


expGenes <- read.csv('../Project/ML/cnv_genes_exp.combine.csv', sep = '\t', header = T)
View(genes) # has genes and do lima analysis

survival <- read.csv('../Project/ML/TCGA-BLCA.survival.tsv', sep = '\t', header = T)
View(survival) # has genes and do lima analysis

survival$sample <- gsub("-", ".", survival$sample)

gene_names <- expGenes$genes

df <- expGenes[1:ncol(expGenes)-1]

df_T <- as.data.frame(t(df))
df_T

colnames(df_T) <- gene_names
View(df_T)

# ----------------------------- Preprocessing ----------------------------
#genes = as.matrix(genes)

#rownames(genes) = genes[,1]
#genes_exp = genes[,2:ncol(genes)]
#genes_dimnames = list(rownames(genes_exp),colnames(genes_exp))
#genes = matrix(as.numeric(as.matrix(genes_exp)), nrow = nrow(genes_exp),dimnames = genes_dimnames)

# -------------------------- merge tables --------------------------------
#mergedTabs <- merge(survival, df_T, by.x = "")
mergedTabs <- cbind(survival, df_T)


# -------------------------- random forest -------------------------------
data.imputed <- rfImpute(sample ~ ., data = survival, iter = 6)


# ----------------------- logistic regression ----------------------------
model <- glm(survival$sample ~ survival$OS, data = survival, family = "binomial")


# --------------------------- other plots --------------------------------
pca_res <- prcomp(survival, scale. = TRUE)

autoplot(df_T_matrix) # this works too
#autoplot(pca_res, data = df_T_matrix, colour = 'Species')

# heatmap
df_T_matrix <- as.matrix(df_T[3:15])
df_T_matrix
heatmap(df_T_matrix)
