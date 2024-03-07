###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

###################################################
### Validation using human bulk RNA-seq dataset ###
###################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

### Load clinical information and bulk RNA-seq profile data
clinical_info <- read.csv('~/R/invivo/external_NM_220623/EGAD00001008128_clinical_data_20210913.csv')
RNA_mtx <- read.table('~/R/invivo/external_NM_220623/EGAD00001008128_RNA_reverse_readcount.txt', header = T, row.names = 1)
### The data was from the article: Molecular correlates of clinical response and resistance to atezolizumab in combination with bevacizumab in advanced hepatocellular carcinoma
### https://pubmed.ncbi.nlm.nih.gov/35739268/
### Dataset: EGAD00001008128[https://ega-archive.org/datasets/EGAD00001008128]

### We only analyzed the pre-treatment data from group F2
clinical_info <- clinical_info[clinical_info$Treatment.Group == 'F2' & clinical_info$Visit == 'Pre-treatment',]
RNA_mtx <- RNA_mtx[,colnames(RNA_mtx) %in% clinical_info$anon_sampleId]

### Reorder clinical data
reorder_index <- match(colnames(NM_RNA), rownames(NM_clinical))
NM_clinical <- NM_clinical[reorder_index,]

###############################
### Step 1, Data filtration ###
###############################

### Filtering the samples
### Leave only CR, PR, SD and PD samples (which contains response information toward ICI therapy)
filter_pat <- clinical_info$Confirmed.Response_IRF %in% c('CR','PR','SD','PD')
### 43 samples in group F2
RNA_mtx <- RNA_mtx[,filter_pat]
clinical_info <- clinical_info[filter_pat,]

### Filtering the genes
### We filtered the rows (features) with available gene symbol
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_IDs <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"), values = gsub('\\..*','',rownames(NM_RNA)), mart = mart)

### Filter the RNA_mtx
filter <- gsub("\\..*", "", rownames(RNA_mtx)) %in% gene_IDs$ensembl_gene_id
RNA_mtx <- RNA_mtx[filter,]
### RNA_mtx [60466, 43]

##################################
### Step 2, Data preprocessing ###
##################################

### We performed count normalization using DEseq2
### Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = RNA_mtx, colData = clinical_info, design = ~ Confirmed.Response_IRF)
### Median of ratios method of normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

### save filtered and normalized data
write.table(normalized_counts, file = paste0(work_path, "Human_bulk_normalized_counts.txt", sep="\t", quote=F, col.names=NA))
            
####################
### Step 3, GSVA ###
####################

### We collected the ENSG number of Runx2 modules
### RUNX2 related genes include IL18RAP, KLRK1, STAT3, RBPJ, LGALS3, RUNX2, NRP1, CTLA4
### ENSG: IL18RAP(ENSG00000115607.9), KLRK1(ENSG00000213809.9), STAT3(ENSG00000168610.15), RBPJ(ENSG00000168214.20), 
### ENSG: LGALS3(ENSG00000131981.16), RUNX2(ENSG00000124813.23), NRP1(ENSG00000099250.18), CTAL4(ENSG00000163599.17)
RUNX2_gene <- c('ENSG00000115607.9', 'ENSG00000213809.9', 'ENSG00000168610.15', 'ENSG00000168214.20', 'ENSG00000131981.16', 'ENSG00000124813.23', 'ENSG00000099250.18', 'ENSG00000163599.17'))

### We used normalized matrix for GSVA
### Then, to mimic T cell expression profile, we normalized the gene expression profile using CD8A expression level
cd8a_factor <- as.data.frame(normalized_counts[gsub('\\..*','',rownames(normalized_counts)) %in% 'ENSG00000153563',])
cd8a_factor <- cd8a_factor[rep(1, 60466),]

### GSVA
gsva.runx2 <- gsva(as.matrix(normalized_counts/cd8a_factor), list(RUNX2_gene), kcdf = 'Poisson', verbose=FALSE)
### Please refer to Section 4 for the main figure code 

######################################
### Step 4, RUNX2 expression level ###
######################################

RUNX2_normalized_counts <- normalized_counts[rownames(normalized_counts) %in% RUNX2_gene,]

### To mimic the T cell expression profile, we simply normalized the RUNX2 module using CD8A expression level
cd8a_factor <- as.data.frame(normalized_counts[gsub('\\..*','',rownames(normalized_counts)) %in% 'ENSG00000153563',])
### Repeat 8 times to fit the Runx2 module matrix
cd8a_factor <- cd8a_factor[rep(1, 8),]
RUNX2_CD8A <- RUNX2_normalized_counts/cd8a_factor

### Data scaling
RUNX2_CD8A.scale <- t(scale(t(RUNX2_CD8A)))

### Heatmap visualization [Supp. Fig. 9b]
### Make colAnn_df
colAnn_df <- as.data.frame(NM_clinical$Confirmed.Response_IRF)
colnames(colAnn_df) <- c('Patient_type')
rownames(colAnn_df) <- rownames(NM_clinical)
for(i in 1:nrow(colAnn_df)){
  if(colAnn_df$Patient_type[i] %in% c('CR','PR')){
    colAnn_df$Response[i] <- 'Response'
  }else{
    colAnn_df$Response[i] <- 'Non-Response'
  }
}

### Reorder patient type
colAnn_df <- colAnn_df[order(match(colAnn_df$Patient_type, c('CR','PR','SD','PD'))),]
reorder_index_pat <- match(rownames(colAnn_df), colnames(RUNX2_CD8A.scale))
RUNX2_CD8A.scale <- RUNX2_CD8A.scale[,reorder_index_pat]
RUNX2_rownames <- c('IL18RAP', 'CTLA4','STAT3','RUNX2','NRP1','KLRK1','LGALS3','RBPJ')
rownames(RUNX2_CD8A.scale) <- RUNX2_rownames

### Color list
colours <- list(
  Patient_type = c('CR' = '#90DB81','PD' = '#ED7268','PR' = '#C5E6C0','SD' = '#F4B9B5'),
  Response = c('Response' = '#90DB81', 'Non-Response' = '#ED7268'))

### Heatmap annotation
colAnn <- HeatmapAnnotation(
  df = colAnn_df,
  which = 'row', ### 'col' (samples) or 'row' (gene) annotation?
  na_col = 'white', ### default colour for any NA values in the annotation data-frame, 'ann'
  col = colours,
  annotation_height = 0.6,
  annotation_width = unit(1, 'cm'),
  gap = unit(1, 'mm'),
  annotation_legend_param = list(
    Patient_type = list(
      at = c('CR', 'PR', 'SD', 'PD'),
      nrow = 4, ### number of rows across which the legend will be arranged
      title = 'Patient type',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')),
    Response = list(
      at = c('Response','Non-Response'),
      nrow = 2,
      title = 'Response',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold'))))

### Plot
Heatmap(t(RUNX2_CD8A.scale),
        cluster_row_slices = FALSE,
        name = 'Z-score',
        col = colorRamp2(c(-2,-1,0,1,2), c('#377EB8',"#BAD6E8", 'white','#FDCEB5','#D4524E')),
                
        ## Row (gene) parameters
        cluster_rows = FALSE,
        show_row_dend = TRUE,
        row_title = '',
        row_title_side = 'left',
        row_title_gp = gpar(fontsize = 16,  fontface = 'bold'),
        row_title_rot = 90,
        show_row_names = F,
        row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
        row_names_side = 'left',
        row_dend_width = unit(25,'mm'),
                
        ### Column (sample) parameters
        cluster_columns = FALSE,
        show_column_dend = TRUE,
        column_title = 'Gene symbol\n',
        column_title_side = 'bottom',
        column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
        column_title_rot = 0,
        show_column_names = T,
        column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
        column_names_rot = 45,
        column_names_max_height = unit(10, 'cm'),
        column_dend_height = unit(25,'mm'),
                
        ### Cluster methods for rows and columns
        clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
        clustering_method_columns = 'ward.D2',
        clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
        clustering_method_rows = 'ward.D2',
                
        ### Specify top and bottom annotations
        right_annotation = colAnn,
        width = 8*unit(5, "mm"), 
        
        ### Parameters for the colour-bar that represents gradient of expression
        heatmap_legend_param = list(
          color_bar = 'continuous',
          legend_direction = 'horizontal',
          legend_position = 'bottom',
          legend_width = unit(5, 'cm'),
          legend_height = unit(5.0, 'cm'),
          title_position = "leftcenter",
          title_gp=gpar(fontsize = 12, fontface = 'bold'),
          labels_gp=gpar(fontsize = 12, fontface = 'bold')))

### save file
write.csv(gsva.runx2, paste0(work_path, "Human_bulk_gsva.csv"))
save(gsva.runx2, colAnn_df, file = past0(work_path, "Human_bulk.RData"))
