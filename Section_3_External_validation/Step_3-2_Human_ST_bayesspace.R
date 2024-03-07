###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##############################################################
### Validation using human spatial transcriptomics dataset ###
##############################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

### Load ST data
load(paste(work_path, "Human_ST_annotated.RData"))

###################################
### Step 1, BayesSpace analysis ###
###################################

### We conducted both the spatial and feature enhancement through BayesSpace
### We only select several genes undergo feature enhancement
### Including the Runx2 module, exhausted module and T cell markers [CD3D, CD3E and CD8A]
feature <- c('RUNX2','CTLA4','RBPJ','IL18RAP','LGALS3','KLRK1','NRP1','STAT3','CD3D','CD3E','CD8A', 'PDCD1', 'TOX','TIGIT','LAG3', 'HAVCR2')

feature_enhancement <- function(data, image_pos, q, feature){
  diet.seurat = Seurat::DietSeurat(data, graphs = "Spatial") ### Slim down Seurat obj prior to conversion
  sce = as.SingleCellExperiment(diet.seurat) ### Convert seurat to SCE
  colData(sce) = cbind(colData(sce), image_pos) ### Add spatial info to SCE
  
  ### BayesSpace Workflow
  sce = spatialPreprocess(sce, platform = "Visium", skip.PCA = F, log.normalize = T) ### Add BayesSpace metadata, without messing with PCA/logcounts
  sce = spatialCluster(sce, nrep = 1000, burn.in = 100, q = q) ### Quickly cluster via BayesSpace
  
  ### Spatial enhancement  
  sce.enhanced <- spatialEnhance(sce, q=q, platform="Visium", d=20,
                                 model="t",
                                 nrep=1000, burn.in=100,
                                 save.chain=TRUE)

  ### Feature enhancement
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce, assay.type = 'logcounts', feature_names = feature, model = c("xgboost"))
  return(sce.enhanced)
}


sce.enhanced1 <- feature_enhancement(P1T, P1T@images[["image_P11_T"]]@coordinates, 14, feature)
sce.enhanced3 <- feature_enhancement(P3T, P3T@images[["image_P15_T"]]@coordinates, 17, feature)
sce.enhanced5 <- feature_enhancement(P5T, P5T@images[["image"]]@coordinates, 18, feature)
sce.enhanced7 <- feature_enhancement(P7T, P7T@images[["image"]]@coordinates, 15, feature)
sce.enhanced8 <- feature_enhancement(P8T, P8T@images[["image"]]@coordinates, 20, feature)
sce.enhanced9 <- feature_enhancement(P9T, P9T@images[["image"]]@coordinates, 22, feature)
sce.enhanced10 <- feature_enhancement(P10T, P10T@images[["image"]]@coordinates, 20, feature)
sce.enhanced11 <- feature_enhancement(P11T, P11T@images[["image"]]@coordinates, 22, feature)
