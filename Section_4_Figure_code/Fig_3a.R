########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 3a ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS("./murine_tcell_Runx28.rds")
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

### Set default assay as TF
DefaultAssay(HCC.tcell) <- "TF"
DefaultAssay(tcell) <- "TF"

### Extract DE-TF
TF.markers.HCC.tcell <- FindMarkers(HCC.tcell, only.pos = TRUE, ident.1 = c('1','4'), ident.2 = c('2','7'), slot = 'count', min.pct = 0.1, logfc.threshold = 0.02)

TF.markers.invivo <- TF.markers.invivo[order(TF.markers.invivo$avg_log2FC, decreasing = T),]
TF.markers.invivo <- TF.markers.invivo[1:10,]
