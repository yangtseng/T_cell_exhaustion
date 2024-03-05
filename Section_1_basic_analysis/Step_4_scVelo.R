#################################
### Section 1, Basic analysis ###
#################################

###############################
### RNA velocity via scVelo ###
###############################

### The scVelo is running under python 
### Here, we converted original Seurat object into h5ad for scanpy and scVelo in python

###################################
### Step 1, write H5Seurat file ###
###################################

SaveH5Seurat(HCC.tcell, '~/Desktop/HCC_tcell.H5Seurat', overwrite = TRUE)
