###################################################################
### Section 3, External validation of Runx2 in exhausted T cell ###
###################################################################

##################################################
### Validation using murine microarray dataset ###
##################################################

set.seed(1234)
work_path = "./external_validation/"
source("requirements.R")

##################################
### Step 1, Data preprocessing ###
##################################

### Load microarray data
murine_data <- ReadAffy(celfile.path = paste0(work_path, "murine_data"))

### Preprocessing of microarray data
eset.rma <- rma(murine_data)
exprs <- exprs(eset.rma)

### Combination with probeid
probeid <- rownames(exprs)
exprs <- cbind(probeid, exprs)

### Annotation
### Of note, here we used htmg430pm.db as annotation database
ID <- AnnotationDbi::select(htmg430pm.db, rownames(exprs), keytype = 'PROBEID', columns = c('SYMBOL','ENTREZID','GENENAME'))

### Extraction of Runx2 and its related genes from microarray data
Runx2_module <- c("Runx2", "Lgals3", "Nrp1", "Ctla4", "Klrk1", "Rbpj", "Stat3", "Il18rap")
ID_runx2 <- ID[ID$SYMBOL %in% Runx2_module,]
eset.rma_runx2 <- eset.rma[rownames(eset.rma@assayData[["exprs"]]) %in% ID_runx2$PROBEID,]

### We summed up the probes annotated to same gene symbol
runx2_exprs <- eset.rma_runx2@assayData[['exprs']]
ID_runx2 <- ID_runx2[order(ID_runx2$SYMBOL),]
runx2_exprs <- runx2_exprs[match(ID_runx2$PROBEID, rownames(runx2_exprs)),]

### Get colSums of each gene to replace duplicated probeIDs
Ctla4 <- runx2_exprs[1,]
Il18rap <- colSums(runx2_exprs[2:3,])
Klrk1 <- runx2_exprs[4,]
Lgals3 <- runx2_exprs[5,]
Nrp1 <- colSums(runx2_exprs[6:9,])
Rbpj <- colSums(runx2_exprs[10:12,])
Runx2 <- colSums(runx2_exprs[13:15,])
Stat3 <- colSums(runx2_exprs[16:18,])

### Generate new matrix of colSummed data
runx2_exprs_sum <- rbind(Ctla4,Il18rap,Klrk1,Lgals3,Nrp1,Rbpj,Runx2,Stat3)

######################################
### Analysis of microarray dataset ###
######################################

# Basic analysis of microarray dataset
pData(eset.rma)
strain <- c("ctrl","ctrl","ctrl","ctrl","treat","treat","treat")
### The microarray data contains four control samples and three treatment samples

design <- model.matrix(~factor(strain))

fit <- lmFit(runx2_exprs_mix, design)
fit <- eBayes(fit)

### options(digits=2)
res <- topTable(fit, number = 'Inf', adjust.method="fdr", coef=2)
res$probe <- rownames(res)

### save file
write.csv(res, paste0(work_path, 'DEG_sum.csv')

### Data scaling
runx2_exprs_sum.scaled = t(scale(t(runx2_exprs_sum)))
colnames(runx2_exprs_sum.scaled) <- c('M5.1\nPlacebo','M12.1\nPlacebo','M14.4\nPlacebo','M17.1\nPlacebo','M12.4\nAnti-PD1','M16.1\nAnti-PD1','M17.4\nAnti-PD1')
### If you are looking for the code of heatmap visualization, please refer to Section_4_Code_for_main_figure

#####################
### GSVA analysis ###
#####################

### We first created gene symbol-wise matrix from pre-processed microarray expression matrix
### ID cleanerance
ID_all <- ID[!is.na(ID$SYMBOL),] ### 40739 obs.
ID_all <- ID_all[!duplicated(ID_all$PROBEID),]
### Reorder by PROBEID
ID_all <- ID_all[order(ID_all$PROBEID),]

### eset.rma cleanerance
eset.rma_all <- eset.rma[rownames(eset.rma@assayData[["exprs"]]) %in% ID_all$PROBEID,]
reorder_index <- match(ID_all$PROBEID, rownames(eset.rma_all@assayData[["exprs"]]))
eset.rma_allmat <- eset.rma_all@assayData[["exprs"]][reorder_index,]

### To deal with the duplicated symbol, we used colSums to get the max value of each symbol and removed the duplicated one
### colSums for the duplicated symbol
dup <- c()
for(i in 1:nrow(ID_all)){
  dup <- c(i)
  for(j in i:nrow(ID_all)){
    if(ID_all$SYMBOL[j] == ID_all$SYMBOL[i]){
      dup <- c(dup, (j))
    }
  }
  if(length(dup) > 1){
    eset.rma_allmat[i,] <- colSums(eset.rma_allmat[dup,])
  }
}

eset.rma_allmat.sum <- eset.rma_allmat

# remove duplicated genes and assign gene symbol to replace probeid
remove_probeid <- function(data, id){
  filter <- duplicated(id$SYMBOL)
  data <- data[!filter,]
  id <- id[!filter,]
  rownames(data) <- id$SYMBOL
  return(data)
}

allmat.sum <- remove_probeid(eset.rma_allmat.sum, ID_all)

### Calcualte the GSVA score of Runx2 module
gsva.runx2_sum <- gsva(allmat.sum, list(Runx2_module), verbose=FALSE)
strain_runx2 <- c('Placebo','Placebo','Placebo','Placebo','Anti-PD1','Anti-PD1','Anti-PD1')
### If you are looking for the code of boxplot visualization, please refer to Section_4_Code_for_main_figure

### save RData for visualization
save(runx2_exprs_sum.scaled, gsva.runx2_sum, Runx2_module, strain_runx2, file = paste0(work_path, "murine_microarray.RData"))
