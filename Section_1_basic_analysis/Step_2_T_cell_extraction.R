#################################
### Section 1, Basic analysis ###
#################################

###########################################################################################################
### T cell extraction, cell composition, T cell subtype annotation and differential expression analysis ###
###########################################################################################################

##########################################################
### T cell extraction from pre-processed seurat object ###
##########################################################
set.seed(1234)
work_path = "./"
source("requirements.R")

load("murine_doubletremoval2.rds")
### It will load a pre-processed seurat object of all cells from section 1, step 1

### Add experimental information based on the murine model design
HCC_time <- as.character(HCC.s@meta.data[["orig.ident"]])
HCC_time <- substring(HCC_time, 1, 4)
for(i in 1:length(HCC_time)){
  if(HCC_time[i] == "P136"){
    HCC_time[i] <- "57_day"
  }else{
    HCC_time[i] <- "20_day"
  }
}

### Add the experimental information to seurat object 
HCC.s <- AddMetaData(HCC.s, metadata = as.factor(HCC_time), col.name = 'time')

### T cell extraction based on the expression level of Cd3e, Cd8a
### Visualization of Cd3e and Cd8a expression level [Supp. Fig. S2]
FeaturePlot(HCC.s, features = c('Cd3e','Cd8a'), ncol = 2)

### Cluster-based extraction of T cells from seurat object
HCC_s.split <- SplitObject(HCC_s, split.by = 'seurat_clusters')
HCC.tcell <- merge(HCC_s.split[["1"]], y = c(HCC_s.split[["2"]], HCC_s.split[["3"]], HCC_s.split[["4"]], HCC_s.split[["7"]], HCC_s.split[['8']]), project = 'HCC.tcell')
### Cluster 1, 2, 3, 4, 7 and 8 are considered as T cell clusters

### Pre-processing of T cell only
HCC.tcell <- NormalizeData(HCC.tcell, verbose = F)
HCC.tcell <- FindVariableFeatures(HCC.tcell, selection.method = 'vst', nfeatures = 3000)
HCC.tcell <- ScaleData(HCC.tcell, verbose = F)
HCC.tcell <- RunPCA(HCC.tcell, npcs = 50, verbose = F)

### We removed the batch effect through Harmony
HCC.tcell <- RunHarmony(HCC.tcell, 'time')

#UMAP and TSNE
HCC.tcell <- RunUMAP(HCC.tcell, reduction = "harmony", dims = 1:20)
HCC.tcell <- RunTSNE(HCC.tcell, reduction = 'harmony', dims = 1:20)
HCC.tcell <- FindNeighbors(HCC.tcell, reduction = "harmony", dims = 1:20)
HCC.tcell <- FindClusters(HCC.tcell, resolution = 0.5, algorithm = 4, method = 'igraph')


