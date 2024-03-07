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
### RNA_mtx [60512, 43]

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
write.table(normalized_counts, file = paste0(work_path, "Human_bulk_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

##############################
### Step 3, Basic analysis ###
##############################

### We collected the ENSG number of Runx2 modules
### RUNX2 related genes include IL18RAP, KLRK1, STAT3, RBPJ, LGALS3, RUNX2, NRP1, CTLA4
### ENSG: IL18RAP(ENSG00000115607), KLRK1(ENSG00000213809), STAT3(ENSG00000168610), RBPJ(ENSG00000168214), 
### ENSG: LGALS3(ENSG00000131981), RUNX2(ENSG00000124813), NRP1(ENSG00000099250), CTAL4(ENSG00000163599)
RUNX2_gene <- c('ENSG00000115607', 'ENSG00000213809', 'ENSG00000168610', 'ENSG00000168214', 'ENSG00000131981', 'ENSG00000124813', 'ENSG00000099250', 'ENSG00000163599')
RUNX2_normalized_counts <- normalized_counts[gsub('\\..*','',rownames(normalized_counts)) %in% RUNX2_gene,]

### To mimic the T cell expression profile, we simply normalized the gene expression profile using CD8A expression level
cd8a_factor <- as.data.frame(normalized_counts[gsub('\\..*','',rownames(normalized_counts)) %in% 'ENSG00000153563',])
### Repeat 8 times to fit the Runx2 module matrix
cd8a_factor <- cd8a_factor[rep(1, 8),]
RUNX2_CD8A <- RUNX2_normalized_counts/cd8a_factor

### Data scaling
RUNX2_CD8A.scale <- t(scale(t(RUNX2_CD8A)))


