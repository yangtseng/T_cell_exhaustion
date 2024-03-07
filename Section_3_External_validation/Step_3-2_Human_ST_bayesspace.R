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

Feature_enhancement <- function(data, image_pos, q, feature){
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

### Enhancement on each sample
sce.enhanced1 <- Feature_enhancement(P1T, P1T@images[["image_P11_T"]]@coordinates, 14, feature)
sce.enhanced3 <- Feature_enhancement(P3T, P3T@images[["image_P15_T"]]@coordinates, 17, feature)
sce.enhanced5 <- Feature_enhancement(P5T, P5T@images[["image"]]@coordinates, 18, feature)
sce.enhanced7 <- Feature_enhancement(P7T, P7T@images[["image"]]@coordinates, 15, feature)
sce.enhanced8 <- Feature_enhancement(P8T, P8T@images[["image"]]@coordinates, 20, feature)
sce.enhanced9 <- Feature_enhancement(P9T, P9T@images[["image"]]@coordinates, 22, feature)
sce.enhanced10 <- Feature_enhancement(P10T, P10T@images[["image"]]@coordinates, 20, feature)
sce.enhanced11 <- Feature_enhancement(P11T, P11T@images[["image"]]@coordinates, 22, feature)

###############################################
### Step 2, Generation of enhanced spot mtx ###
###############################################

### We combined each subspot and generate enhanced spot matrix to put back to seurat object for further analysis
EnhSpot <- function(seurat_obj, sce.enhanced){
  data <- sce.enhanced@assays@data@listData[["logcounts"]]
  data <- data[!rowSums(is.na(data)),]
  max <- ncol(data)/6
  spot.mtx <- data.frame(matrix(NA, nrow = nrow(data), ncol = max))
  for(j in 1:nrow(data)){
    for(i in 1:max){
      spot.mtx[j,i] <- sum(data[j,i], data[j,i+max],data[j,i+2*max], data[j,i+3*max],data[j,i+4*max], data[j,i+5*max]) 
    }
  }

  ### Assign the rownames
  rownames(spot.mtx) <- rownames(data)
  ### Order with orginal cell name
  colnames(spot.mtx) <- seurat_obj@assays[["Spatial"]]@data@Dimnames[[2]][order(seurat_obj@assays[["Spatial"]]@data@Dimnames[[2]])]
  
  # Assign back to seurat object
  seurat_obj[['Enhanced']] <- CreateAssayObject(counts = spot.mtx)
  return(seurat_obj)
}

P1T <- EnhSpot(P1T, sce.enhanced1)
P3T <- EnhSpot(P3T, sce.enhanced3)
P5T <- EnhSpot(P5T, sce.enhanced5)
P7T <- EnhSpot(P7T, sce.enhanced7)
P8T <- EnhSpot(P8T, sce.enhanced8)
P9T <- EnhSpot(P9T, sce.enhanced9)
P10T <- EnhSpot(P10T, sce.enhanced10)
P11T <- EnhSpot(P11T, sce.enhanced11)

### save files
save(P1T, P3T, P5T, P7T, P8T, P9T, P10T, P11T, file = paste0(work_path, "Human_ST_bayesspace.RData"))
