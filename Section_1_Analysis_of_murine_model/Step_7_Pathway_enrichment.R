###########################################
### Section 1, Anlaysis of murine model ###
###########################################

###################################
### Pathway enrichment analysis ###
###################################

set.seed(1234)
work_path = "./"
source("requirements.R")

HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_TF6.rds"))

###################################################################
### Step 1, AUC score calculation of Hallmark and KEGG gene set ###
###################################################################

### We used AUCell to calculate the enrichment score of pathways for each cell
### We first built the rank via the gene expression profile
CellRank <- AUCell_buildRankings(as.matrix(HCC.tcell@assays$RNA@data))

### Get mouse gene sets data
### Here, we used KEGG and Hallmark gene sets
gs.go.h <- getGeneSets(species = "Mus musculus", library = 'H')
gs.go.kegg <- getGeneSets(species = 'Mus musculus', library = 'C2', subcategory = 'KEGG')

### Calculate the gene set AUC score for Hallmark gene set
cells_AUC_h <- AUCell_calcAUC(gs.go.h, CellRank, aucMaxRank=nrow(CellRank)*0.05)
save(cells_AUC_h, file = paste0(work_path, 'cells_AUC_h.RData'))

### Calculate the gene set AUC score for KEGG gene set
cells_AUC_kegg <- AUCell_calcAUC(gs.go.kegg, CellRank, aucMaxRank=nrow(CellRank)*0.05)
save(cells_AUC_kegg, file = paste0(work_path, 'cells_AUC_kegg.RData'))

################################################
### Step 2, Exploring the enrichment results ###
################################################

### We first added the enrichment results to seurat object
HCC.tcell[['KEGG']] <- CreateAssayObject(cells_AUC_kegg@assays@data@listData[["AUC"]])
HCC.tcell[['Hallmark']] <- CreateAssayObject(cells_AUC_h@assays@data@listData[["AUC"]])

### Next, we identified the activated pathways in specific exhausted T cell clusters (cluster 1 and 4) compared to effector T cells
### On the other hand, we also revealed the pathways activated in specific exhausted T cell clusters (cluster 1 and 4) compared to the rest of T cells
### Please refer to section 4 for the visualization code of main figures
PE1 <- FindMarkers(HCC.tcell, assay = 'KEGG', only.pos = T, logfc.threshold = 0.01, test.use = 't', ident.1 = c('1','4'), ident.2 = c('2', '7'))
PE2 <- FindMarkers(HCC.tcell, assay = 'KEGG', only.pos = T, logfc.threshold = 0.01, test.use = 't', ident.1 = c('1','4'))

### To explore the fold enrichment, we calculated the average AUC score for each cluster
### We then calculated the fold enrichment value by simply divide the sum of cluster 1 and 4 to cluster 2 and 7
AVG <- AverageExpression(HCC.tcell, assays = "KEGG")[["KEGG"]]
AVG$FE <- (AVG$`1` + AVG$`4`)/(AVG$`2`+ AVG$`7`)

### We merged PE1 and PE2 to extract the pathways that is specific enriched in cluster 1 and 4 
PE <- PE1[rownames(PE1) %in% rownames(PE2),]
AVG <- AVG[rownames(AVG) %in% rownames(PE),]

PE <- merge(PE, AVG, by = 'row.names')
 
### save files
saveRDS(HCC.tcell, paste0(work_path, "murine_tcell_pathway7.rds"))
write.csv(PE, paste0(work_path, "Pathway_enrichment_Tex.csv"))
