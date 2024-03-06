###########################################
### Section 1, Anlaysis of murine model ###
###########################################

###########################################################################
### We revealed relationship between Runx2 and pseudotime in this script###
###########################################################################

### Packages
import scanpy as sc
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

### Basic settings
work_path = "./"

### Load H5ad file
HCC_T_scvelo = sc.read_h5ad("HCC_tcell_scvelo.h5ad")

### Extract gene expression matrix 
df = sc.get.obs_df(HCC_T_scvelo, ['Runx2', 'Nrp1','Stat3','Ctla4','Il18rap','Rbpj','Klrk1','Lgals3', "seurat_clusters", "latent_time"])

### Group cells into cluster 1 and 4 (c1_4) or rest of T cells (rest_of_T)
df['condition'] = np.where((df["seurat_clusters"] == 0) | (df['seurat_clusters'] == 3), 'c1_4','rest_of_T')

### Stack the dataframe by melt()
genes = ['Runx2','Ctla4','Il18rap','Klrk1','Lgals3','Nrp1','Rbpj','Stat3']
ids = ['seurat_clusters','latent_time','condition']

df_s = pd.melt(df, id_vars=ids, value_vars=genes)

### Visualization [Fig. 3a]
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", palette=('#9DD1C7','#E58C82'), rc=custom_params, font_scale = 3)
g = sns.lmplot(data= df_s, x="latent_time", y= 'value', hue="condition", col = 'variable', col_wrap=4,
               sharey= False, order=3, scatter=False, aspect=1.2, legend=None, height=6, line_kws={'linewidth':4}).set(xlabel = '', ylabel = '')
g.fig.supxlabel('Latent time')
g.fig.supylabel('Relative expression value')

### SAve the plot as PNG file
plt.savefig("Runx2_pseudotime.png")