###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################
### Transcription factor prediction ###
#######################################

################
### GRNboost ###
################

### packages
import pandas as pd

from distributed import Client, LocalCluster
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

in_file  = 'int/1.1_exprMatrix_filtered_t.txt'
tf_file  = 'int/1.1_inputTFs.txt'
out_file = 'int/1.1_grn_output.tsv'

### ex_matrix wass a DataFrame with gene names as column names
ex_matrix = pd.read_csv(in_file, sep='\t')

### tf_names was read using a utility function included in Arboreto
tf_names = load_tf_names(tf_file)

### We instantiated a custom Dask distributed Client
client = Client(LocalCluster())

### We computed the GRN
network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names,
                    client_or_address=client)

### We saved the GRN to file
network.to_csv(out_file, sep='\t', index=False, header=False)

