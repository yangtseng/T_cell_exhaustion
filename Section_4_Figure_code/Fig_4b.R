########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

### Set default assay to chromVAR
DefaultAssay(tcell) <- 'chromvar'

### Save

FeaturePlot(
  object = tcell,
  features = "MA0511.2", ### Runx2
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  order = T
)
