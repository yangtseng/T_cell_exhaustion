#################################
### Section 1, Basic analyssi ###
#################################

###############################################
### Data preprocessing and doublets removal ###
###############################################

### This is the section of murine model analysis which can be downloaded in raw and processed forms from the NCBI Gene Expression Omnibus under accession number [NUMBER]
### If you are looking for the cell-line like model analysis, please refer to the section 3

###############################################################
### Step 1, read the raw files and remove low quality cells ###
###############################################################
set.seed(1234)
work_path = "./"

load("murine_raw.RData")
All <- cbind(P136A1, P136A2, P136B1, P136B2, P138A1, P138A2, P138B1, P138B2)

### We performed a three-steps quality control using scater for all cells 
All.qc <- perCellQCMetrics(All, subsets = list(Mito=grep("mt-", rownames(All))))

### First, we removed cells which are outlier (3 nmad) and mito_percent > 5%
outlier <- isOutlier(metric = All.qc$detected, nmad = 3, log = T)
mito_filter <- All.qc$subsets_Mito_percent >= 5

filter <- outlier | mito_filter
All <- All[,!filter]

### Second, we removed low expressed genes and cells
gene_filter <- rowSums(All > 0) >= 10
cell_filter <- colSums(All > 0) >= 100

All <- All[gene_filter, cell_filter]

### Third, we averaged the duplicated gene_symbol and remove the rest of them
gene_symbol <- All@Dimnames[[1]]
gene_symbol <- gene_symbol[!duplicated(gene_symbol)]

All_1 <- data.table(as.matrix(All))
All_1 <- cbind(All@Dimnames[[1]], All_1)

All_1 <- All_1[,lapply(.SD, mean),by = V1]
All_2 <- All_1[,2:ncol(All_1)]

rownames(All_2) <- gene_symbol

###########################################################################
### Step 2, preprocess the single-cell RNA-seq dataset with Seurat in R ###
###########################################################################

### We first created the Seurat object based on filtered data
HCC <- CreateSeuratObject(counts = All_2, project = "HCC_scRNA")
HCC <- NormalizeData(HCC, verbose = FALSE)
HCC <- FindVariableFeatures(HCC, selection.method = "vst", nfeatures = 3000)

### Data scaling and principle component analysis
HCC <- ScaleData(HCC, verbose = FALSE)
HCC <- RunPCA(HCC, npcs = 50, verbose = FALSE)

### To determine the best elbow point, we calculated the cumulative percentage for each PC 
elbow.plot.info <- ElbowPlot(HCC)
ElbowPlot(HCC, ndims = 50)

pct <- HCC[["pca"]]@stdev / sum(HCC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

### We determined the elbow point based on two criterias
### [1] PC exhibited cumulative percent greater than 90% and % variation associated with the PC as less than 5
criteria1 <- which(cumu > 90 & pct < 5)[1]
### [2] the difference between variation of PC and subsequent PC
criteria2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

### last point where change of % of variation is more than 0.1%.
pcs <- min(criteria1, criteria2)
print(pcs)

### Lastly, we conducted two-dimensional reduction (UMAP/t-SNE) and clustering using UMAP and Leiden algorithm (igraph)
HCC <- RunUMAP(HCC, reduction = "pca", dims = 1:pcs)
HCC <- RunTSNE(HCC, reduction = 'pca', dims = 1:pcs)
HCC <- FindNeighbors(HCC, reduction = "pca", dims = 1:pcs)
HCC <- FindClusters(HCC, resolution = 0.5, algorithm = 4, method = 'igraph')

### Save RDS
saveRDS(HCC, paste0(work_path, "murine_preprocessed1.rds"))

################################
### Step 3, Doublets removal ###
################################

### Split seurat object into 8 subsets by sample
sc.list <- SplitObject(HCC, split.by = 'orig.ident')

### We performed doublets removal using Doubletfinder in R
phe_lt <- lapply(names(sc.list), function(x){
  sc.filt <- sc.list[[x]]
  sc.filt <- FindVariableFeatures(sc.filt)
  sc.filt <- ScaleData(sc.filt)
  sc.filt <- RunPCA(sc.filt, npcs = 50)
  sc.filt <- RunUMAP(sc.filt, dims = 1:pcs)
  
  ### Doubletfinder
  ###########################################
  ### pK Identification (no ground-truth) ###
  sweep.res.sc.filt <- paramSweep_v3(sc.filt, PCs = 1:pcs, sct = FALSE)
  sweep.stats.sc.filt <- summarizeSweep(sweep.res.sc.filt, GT = FALSE)
  bcmvn.sc.filt <- find.pK(sweep.stats.sc.filt)
  pK <- bcmvn.sc.filt$pK[which(bcmvn.sc.filt$BCmetric == max(bcmvn.sc.filt$BCmetric))]
  pK <- as.numeric(as.character(pK))
  ### Homotypic Doublet Proportion Estimate ###
  annotations <- sc.filt@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.05*ncol(sc.filt))  ## Assuming 5% doublet formation rate - tailor for your dataset
  
  ### Run DoubletFinder with varying classification stringencies ###
  sc.filt <- doubletFinder_v3(sc.filt, PCs = 1:pcs, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  ##########################################
})

### Results of Doubletfinder 
kp.1 <- phe_lt[[1]]@meta.data[["DF.classifications_0.25_0.19_577"]]
kp.2 <- phe_lt[[2]]@meta.data[["DF.classifications_0.25_0.03_454"]]
kp.3 <- phe_lt[[3]]@meta.data[["DF.classifications_0.25_0.23_650"]]
kp.4 <- phe_lt[[4]]@meta.data[["DF.classifications_0.25_0.1_701"]]
kp.5 <- phe_lt[[5]]@meta.data[["DF.classifications_0.25_0.02_246"]]
kp.6 <- phe_lt[[6]]@meta.data[["DF.classifications_0.25_0.26_244"]]
kp.7 <- phe_lt[[7]]@meta.data[["DF.classifications_0.25_0.3_428"]]
kp.8 <- phe_lt[[8]]@meta.data[["DF.classifications_0.25_0.01_308"]]

kp <- c(kp.1, kp.2, kp.3, kp.4, kp.5, kp.6, kp.7, kp.8)
kpc <- as.factor(kp)

### Add doubletfinder result to seurat project
HCC <- AddMetaData(HCC, metadata = kpc, col.name = 'Doublet')

### Doublet removal
HCC.split <- SplitObject(HCC, split.by = 'Doublet')
HCC_s <- HCC.split[['Singlet']]

### Save RDS
saveRDS(HCC, paste0(work_path, "murine_doubletremoval2.rds"))

