########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 3c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
gs.go.kegg <- escape::getGeneSets(species = 'Mus musculus', library = 'C2', subcategory = 'KEGG')

### Filter the KEGG pathway based on the result of pathway enrichment 
KEGG <- gs.go.kegg[names(gs.go.kegg) %in% gsub("-", "_", rownames(KEGG_AVG27))]
