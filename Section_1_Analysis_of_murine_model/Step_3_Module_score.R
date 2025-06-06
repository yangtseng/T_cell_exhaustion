###########################################
### Section 1, Anlaysis of murine model ###
###########################################

##########################################################
### Calculation of module score and cell cycle scoring ###
##########################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path,"murine_tcell3.rds"))
### It will load a pre-processed seurat object of T cells from section 1, step 2

################################
### Module score calculation ###
################################

### We calcuated 3 different module scores including effectory, exhaustion and immune checkpoint modules
### The gene sets were collected from other related articles
### We calculated the module score of each cell via AddModuleScore() in seurat

### T cell effectory gene set
eff <- list(c('Ifng','Ccl3','Ccl4','Prf1','Nkg7','Gzmb','Gzmk'))
HCC.tcell <- AddModuleScore(HCC.tcell, features = eff, name = 'effectory')

### T cell exhaustion gene set
exhaust <- list(c('Pdcd1','Cd244a', 'Lag3','Tigit','Eomes','Tox'))
HCC.tcell <- AddModuleScore(HCC.tcell, features = exhaust, name = 'exhaustion')

### T cell checkpoint gene set
cp <-  list(c('Havcr2', 'Tigit', 'Ctla4', 'Pdcd1', 'Lag3', 'Entpd1'))
HCC.tcell <- AddModuleScore(HCC.tcell, features = cp, name = 'checkpoint')

####################################
### Cell cycle score calculation ###
####################################

### We calcuated the cell cycle score following the seurat tutorial
### The gene sets were collected from Tirosh et al, 2015 which could be loaded with Seurat directly
### To fit the murine model, we first converted the gene symbol from human to mice
### Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)  
  humanx <- unique(genesV2[, 2])

  return(humanx)
}

### Cell cycle gene lists for mice
ms.s.genes <- convertHumanGeneList(cc.genes$s.genes)
ms.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)

### Cell cycle score calculation via CellCycleScoring()
HCC.tcell <- CellCycleScoring(HCC.tcell, s.features = ms.s.genes, g2m.features = ms.g2m.genes, set.ident = TRUE)

### save RDS for pre-processed T cell
saveRDS(HCC.tcell, paste0(work_path, "murine_tcell_modulescore4.rds"))
### If you would like to check out the code for the visualization of main figures, please refer to section 4.
