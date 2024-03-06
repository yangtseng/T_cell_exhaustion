###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

######################
### Motif analysis ###
######################

set.seed(1234)
work_path = "./Cellline/"
source("requirements.R")

tcell <- readRDS(paste0(work_path, "Cellline_TF3.rds"))

##########################
### Step 1, Link peaks ###
##########################

### This script was about using Signac to analyze the scATAC-seq data from cell line like model
### We linked the peak and gene via Signac function LinkPeaks

### Set default assay as peaks
DefaultAssay(tcell) <- 'peaks'

### We first computed the GC content for each peak
tcell <- RegionStats(tcell, genome, BSgenome.Mmusculus.UCSC.mm10)

### We linked peaks to genes
tcell <- LinkPeaks(tcell, peak.assay = "peaks", expression.assay = "SCT")

####################################################################
### Step 2, Identification of differentail accessibility regions ###
####################################################################

### Matrix sets
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", all_version = FALSE))

### Change sequence info of the seurat object
tcell@assays[["peaks"]]@ranges@seqinfo@genome[] <- "mm10"

### Add motifs
tcell <- AddMotifs(tcell, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

### Exploring the DARs in cell line like model
DARs <- FindMarkers(tcell, ident.1 = c('Exhausted_T_72hr','Exhausted_T_96hr'), ident.2 = 'Active_T', 
                    only.pos = TRUE, test.use = "LR", min.pct = 0.05, latent.vars = "nCount_peaks")
DARs$name <- rownames(DARs)

### save files 
write.csv(DARs, paste0(work_path, "DARs.csv"))
saveRDS(tcell, paste0(work_path, "Cellline_motif4.rds"))
