###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

########################################################################
### Module score, gene activity and differential expression analysis ###
########################################################################

set.seed(1234)

work_path = "./Cellline/"
source(paste0(work_path, "requirements.R"))

### Here, we loaded the integrated multiome seurat object
tcell <- readRDS("Cellline_multiome1.rds")

##########################################
### Step 1, Calculate the module score ###
##########################################

### Set default assay as SCT
DefaultAssay(tcell) <- "SCT"

### Calculate the module score of exhausted and effectory gene set
### T cell effectory gene set
eff <- list(c('Ifng','Ccl3','Ccl4','Prf1','Nkg7','Gzmb','Gzmk'))
tcell <- AddModuleScore(tcell, features = eff, name = 'effectory')

### T cell exhaustion gene set
exhaust <- list(c('Pdcd1','Cd244a', 'Lag3','Tigit','Eomes','Tox'))
tcell <- AddModuleScore(tcell, features = exhaust, name = 'exhaustion')

###########################################
### Step 2, Calculate the gene activity ###
###########################################

### Set default assay as peaks 
DefaultAssay(tcell) <- "peaks"

### Gene activity calculation
gene.activities <- GeneActivity(tcell)

### Create a new assay for gene activity storage
tcell[["GA"]] <- CreateAssayObject(counts = gene.activities)

### Normalization and scaling
DefaultAssay(exvivo240202) <- "GA"
tcell <- NormalizeData(tcell)
tcell <- ScaleData(tcell)

################################################
### Step 3, Differential expression analysis ###
################################################

### Identifying the differentially expressed genes (DEGs)
tcell_DEGs <- FindAllMarkers(tcell, assay = 'SCT', logfc.threshold = 0.4, only.pos = T)
tcell_DEGs <- tcell_DEGs[unique(tcell_DEGs$gene),]
tcell_DEGs$cluster <- factor(tcell_DEGs$cluster, levels = c('Exhausted_T_96hr',"Exhausted_T_72hr", "Active_T","Naive_T"))
tcell_DEGs <- tcell_DEGs[order(tcell_DEGs$cluster, tcell_DEGs$avg_log2FC, decreasing = T),]

### Identifying the differentially activated genes (based on gene activity)
tcell_DAGs <- FindAllMarkers(tcell, assay = 'GA', logfc.threshold = 0.3, only.pos = T)
tcell_DAGs <- tcell_DAGs[unique(tcell_DAGs$gene),]
tcell_DAGs$cluster <- factor(tcell_DAGs$cluster, levels = c('Exhausted_T_96hr',"Exhausted_T_72hr", "Active_T","Naive_T"))
tcell_DAGs <- tcell_DAGs[order(tcell_DAGs$cluster, tcell_DAGs$avg_log2FC, decreasing = T),]

### Preparing dataframe for main figure [Fig. 2b]
counts <- GetAssayData(tcell, assay="SCT", slot="data")
counts <- as.matrix(counts[rownames(counts) %in% genes, ])

GA_counts <- GetAssayData(tcell, assay="gene_activity", slot="data")
GA_counts <- as.matrix(GA_counts[rownames(GA_counts) %in% tcell_DAGs$gene, ])

### Scaling the matrix
counts_scaled = t(scale(t(counts)))
GA_counts_scaled = t(scale(t(GA_counts)))

### save files
write.csv(tcell_DEGs, paste0(work_path, "cell_line_like_DEGs.csv"))
write.csv(tcell_DAGs, paste0(work_path, "cell_line_like_DAGs.csv"))

saveRDS(tcell, paste0(work_path, "Cellline_multiome2.rds")
save(counts_scaled, GA_counts_scaled, file = paste0(work_path, "Cellline_DEGs.RData"))
