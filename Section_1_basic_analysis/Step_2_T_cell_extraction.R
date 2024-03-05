#################################
### Section 1, Basic analysis ###
#################################

###########################################################################################################
### T cell extraction, cell composition, T cell subtype annotation and differential expression analysis ###
###########################################################################################################

############################################################################
### T cell extraction and visualization from pre-processed seurat object ###
############################################################################
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

### We conducted two-dimensional reduction (UMAP/t-SNE) and clustering using UMAP and Leiden algorithm (igraph)
### We used only 20 PCs for T cell under harmony reduction
HCC.tcell <- RunUMAP(HCC.tcell, reduction = "harmony", dims = 1:20)
HCC.tcell <- RunTSNE(HCC.tcell, reduction = 'harmony', dims = 1:20)
HCC.tcell <- FindNeighbors(HCC.tcell, reduction = "harmony", dims = 1:20)
HCC.tcell <- FindClusters(HCC.tcell, resolution = 0.5, algorithm = 4, method = 'igraph')
### T cells were clustering into 9 clusters

######################################################
### T cell subtype annotation and cell composition ###
######################################################

### T cell subtype annotation
### We used several well-known T cell subtype marker and conducted cluster-based annotation

### Effector T cell: Cd7, Gzmb, Gzmk
### Exhausted T cell: Ctla4, Lag3, Pdcd1
### High IFN response T cell: Ifit1, Ifit3, Isg15
### Memory T cell: Il7r, Tcf7, S1pr1
### Proliferative T cell: Ccna2, Ccnb2, Cdk1

### all fifteen genes were plotted in same figure [Supp. Fig. 3b]
feature <- c('Cd7','Gzmb','Gzmk','Ctla4','Lag3','Pdcd1','Ifit1','Ifit3','Isg15','Il7r','Tcf7', 'S1pr1', 'Ccna2','Ccnb2','Cdk1')
suppfig3b <- FeaturePlot(HCC.tcell, features = feature, ncol = 3, order = F, combine = T, cols = c('grey80',"#D4524E")) & NoAxes() & NoLegend()
suppfig3b <- suppfig3b  & theme(plot.title = element_text(size = 20, face = "italic"))

### T cells were then annotated into 5 distinct subtypes
celltype <- c('Exhausted T cell','Effector T cell', 'Memory T cell', 'Exhausted T cell', 'Exhausted T cell', 'Exhausted T cell', 'Effector T cell', 'Proliferative T cell', 'High IFN response T cell')
names(celltype) <- levels(HCC.tcell)
HCC.tcell <- RenameIdents(HCC.tcell, celltype)
HCC.tcell$cell_type <- HCC.tcell@active.ident

### Visualization [Supp. Fig. 3c]
DimPlot(HCC.tcell, reduction = 'umap', group.by = 'cell_type', cols = c('#f3877f','#c2e8bc','#8ed3c7','#e2e1c3', '#bd89bd'))

### T cell composition
### We revealed the cell composition in each time point
celltype <- as.data.frame(prop.table(table(HCC.tcell$cell_type)))
celltype$Freq <- par(celltype$Freq)
celltype$x <- 1
celltype <- as.data.frame(celltype[order(celltype$Freq),])

### Visualization [Supp. Fig. 3d]
plot <- ggplot(celltype, aes(x = x, y = Freq, fill = Var1)) + geom_bar(stat="identity")

########################################
### Differential expression analysis ###
########################################

### We simply used Seurat FindAllMarkers() for differnetail expression analysis
tcell.marker <- FindAllMarkers(HCC.tcell, logfc.threshold = 0.5, only.pos = T)
write.csv(tcell.marker, paste0(work_path, "tcell_markers.csv"))

### save RDS for pre-processed T cell
saveRDS(HCC.tcell, paste0(work_path, "murine_tcell3.rds"))
