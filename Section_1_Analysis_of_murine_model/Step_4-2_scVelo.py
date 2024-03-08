###########################################
### Section 1, Anlaysis of murine model ###
###########################################

###########################################
### We performed scVelo via this script ###
###########################################

### Packages
import scanpy as sc
import numpy as np
import velocyto as vc
import loompy as lp
import scvelo as scv
import matplotlib.pyplot as plt
import time

### Basic setting
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")

start_time = time.time()

### Load H5ad file
adata = scv.read('HCC_tcell.h5ad')

############################################################
### Step 1, convert the information to fit the loom file ###
############################################################

### Coversion of cluster information 
dict = {0:1, 1:2, 2:3, 3:4, 4:5, 5:6, 6:7, 7:8, 8:9}
adata.obs['clusters'] = adata.obs['seurat_clusters'].values
adata.obs = adata.obs.replace({"clusters": dict})

### Conversion of experimental information
dict2 = {0:"HCC_008_iso", 1:"HCC_008_anti", 2:"HCC_004_iso", 3:"HCC_004_anti"}
adata.obs['name'] = adata.obs['group'].values
adata.obs = adata.obs.replace({'name': dict2})

#################################################################
### Step 2, Data pre-processing with loom file of each sample ###
#################################################################

### Load the loom file including the amount of splicing and non-splicing RNA 
HCC_loom = ["HCC_004_iso.loom", "HCC_004_anti.loom", "HCC_008_iso.loom", "HCC_008_anti.loom"]

### Data cleaning
scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(HCC_loom)

### Merge the pre-processed (in R) h5ad file and loom file
HCC_T_cell = scv.utils.merge(adata, HCC_loom)

### Data Pre-processing
scv.pp.filter_and_normalize(HCC_T_cell)
scv.pp.moments(HCC_T_cell)

########################################
### Step 3, RNA velocity calculation ###
########################################

### Dynamics recovery
scv.tl.recover_dynamics(HCC_T_cell, n_jobs=10)

### Velocity calcuation with dynamical mode
scv.tl.velocity(HCC_T_cell, mode='dynamical')

### Graph construction
scv.tl.velocity_graph(HCC_T_cell, n_jobs = 4)

### Pseudotime estimation [Supp. Fig. 4c]
scv.tl.latent_time(HCC_T_cell)

#############################################
### Step 4, Visualization of RNA velocity ###
#############################################

### Main figure [Fig. 1g]
scvelo = scv.pl.velocity_embedding_stream(HCC_T_cell, basis='umap', size = 50, arrow_size = 1.5, color = 'clusters', figsize = (9,8), save = './main_figure/Fig_1g.png', palette = ['#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd'])

### Plotting for proportions [Supp. Fig. 4a]
scv.pl.proportions(HCC_T_cell, save = 'suppfig4a.png')

### Visualization [Supp. Fig. 4b]
scv.pl.velocity_embedding_grid(HCC_T_cell, basis='umap', palette = ['#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd'])

### Pseudotime [Supp. Fig. 4c]
scv.pl.scatter(HCC_T_cell, color='latent_time', color_map='gnuplot', size=80)

### Save H5ad
HCC_T_cell.write_h5ad(filename = "HCC_tcell_scvelo.h5ad")
