###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

#############################################
### Data integration of 10x multiome data ###
#############################################

set.seed(1234)

work_path = "./"
source(paste0(work_path, "requirements.R"))

### Here, we integrated the single-cell multiome datasets from 4 independent samples
naive_t <- Load(paste0(work_path, "naive_t.rds"))
active_t <- Load(paste0(work_path, "active_t.rds"))
exh_t72 <- Load(paste0(work_path, "exhaust_t72.rds"))
exh_t96 <- Load(paste0(work_path, "exhaust_t96.rds"))

### Merge samples and annotated batches
tcell <- merge(
  naive_t,
  y = c(active_t, exh_t72, exh_t96),
  add.cell.ids = c("naive", "active", "exh72", "exh96"),
  project = "cell_line_like"
)

### Of note, there were totally 41594 cells passed the quality control
###### Naive_T: 2078 cells
###### Active_T: 6873 cells
###### Exhausted_T_72hr: 18348 cells
###### Exhausted_T_96hr: 14295 cells

############################################
### Data preprocessing of scRNA-seq data ###
############################################

### Set default assay to RNA
DefaultAssay(tcell) <- "RNA"

### We normalized the gene expression through sctransform
tcell <- SCTransform(tcell)  

### Dimensional reduction based on scRNA-seq data only
tcell <- RunPCA(tcell)

#############################################
### Data preprocessing of scATAC-seq data ###
#############################################

### Set default assay to peaks
DefaultAssay(tcell) <- "peaks"

### Peaks data pre-processing
tcell <- FindTopFeatures(tcell, min.cutoff = 5)
tcell <- RunTFIDF(tcell)
tcell <- RunSVD(tcell)

#######################################################
### Graph construction of single-cell multiome data ###
#######################################################

### We built a joint neighbor graph using both assays
tcell <- FindMultiModalNeighbors(
  object = tcell,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

### Lastly, we built a joint UMA visualization
tcell <- RunUMAP(
  object = tcell,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

### save the integrated data
saveRDS(tcell, paste0(work_path, "Cellline_multiome6.rds"))
