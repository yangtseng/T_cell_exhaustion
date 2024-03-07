###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

###################################################
### Validation using human bulk RNA-seq dataset ###
###################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

##################################
### Step 1, Data preprocessing ###
##################################

### Load microarray data
murine_data <- ReadAffy(celfile.path = paste0(work_path, "murine_data"))
### The data was from the article: Cabozantinib enhances anti-PD1 activity and elicits a neutrophil-based immune response in hepatocellular carcinoma
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9167725/
### GEO accession number: GSE174770 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174770]

