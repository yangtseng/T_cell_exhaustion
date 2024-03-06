###########################################
### Section 1, Anlaysis of murine model ###
###########################################

############################################################
### Runx2 relationship with markers of T cell exhaustion ###
############################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

### Load seurat object
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Pathway7.rds"))

##########################
### Runx2 module score ###
##########################

### Runx2 module was defined as the Runx2 and its downstream target genes
### Runx2, Ctla4, Nrp1, Lgals3, Rbpj, Stat3, Klrk1 and Il18rap
Runx2_module <- c('Runx2','Ctla4','Nrp1','Lgals3','Rbpj','Stat3','Klrk1','Il18rap')

### Module score calcuation
HCC.tcell <- AddModuleScore(HCC.tcell, features = list(Runx2_sig), name = "Runx2_module")

##########################################################################
### Explore the relationship between Runx2 and exhaustion module score ###
##########################################################################

### Extract the module score from seurat object
Tex_runx2 <- as.data.frame(HCC.tcell@meta.data[["exhaustion"]])
Tex_runx2$runx2 <- HCC.tcell@meta.data[["Runx2_module"]]

### Join with the clsuter information
Tex_runx2$cluster <- HCC.tcell@meta.data[["seurat_clusters"]]

### Join with the T cell subtype information
Tex_runx2$subtype <- HCC.tcell@active.ident
                           
