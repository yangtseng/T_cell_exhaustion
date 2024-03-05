############################################
### Data preprocessing of scRNA-seq data ###
############################################

### Set default assay to RNA
DefaultAssay(tcell) <- "RNA"

### We normalized the gene expression through sctransform
tcell <- SCTransform(tcell)  

### Dimensional reduction based on scRNA-seq data only
tcell <- RunPCA(tcell)

### Lastly, we built a joint UMA visualization
tcell <- RunUMAP(
  object = tcell,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)
