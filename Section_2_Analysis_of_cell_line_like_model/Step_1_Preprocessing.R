###################################################
### Section 2, Anlaysis of cell line like model ###
###################################################

###############################################
### Data preprocessing of 10x multiome data ###
###############################################

### According to the model design, we sequenced the naive, active, exhausted_72h and exhausted_96h T cells separately
### Here, we pre-processed each sample and create seurat object including both scRNA-seq and scATAC-seq data

set.seed(1234)

work_path = "./"
source(paste0(work_path, "requirements.R"))

### We self-defined the preprocessing function and perform on each sample to generate the single-cell multiome seurat object
### The input file includes RNA expression and ATAC fragments
Multiome_preprocessing <- function(sample_path, project){

  ### Load 10x multiome data
  ### scRNA-seq
  counts <- Read10X(data.dir = paste0(sample_path, "filtered_feature_bc_matrix/"))

  ### scATAC-seq
  frag_path <- paste0(sample_path, "atac_fragments.tsv.gz")

  ### Load annotation file
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotation) <- "UCSC"
  ### Here, we used UCSC for annotation

  ### We first created seurat object based on gene expression profile
  tcell <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    project = project
  )

  ### We then added ATAC assay to the seurat object
  tcell[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  
  #############################################
  ### Data preprocessing of scATAC-seq data ###
  #############################################

  ### Set default assay to ATAC
  DefaultAssay(tcell) <- "ATAC"

  ### Calculation of nucleosome signal and TSS enrichment
  tcell <- NucleosomeSignal(tcell)
  tcell <- TSSEnrichment(tcell)
  
  ### Filter out low quality cells
  tcell <- subset(
    x = tcell,
    subset = nCount_RNA > 1000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_ATAC < 100000 &
    nucleosome_signal < 1.1 &
    TSS.enrichment > 5
  )

  ####################
  ### Peak calling ###
  ####################

  ### We called peaks using MACS2 which should be installed in standalone env
  peaks <- CallPeaks(tcell, macs2.path = paste0(work_path, ".local/bin/macs2"))

  ### We removed peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

  ### Furthermore, we quantified counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(tcell),
    features = peaks,
    cells = colnames(tcell)
  )

  ### We created a new assay using the MACS2 peak set and add it to the Seurat object
  ### The peaks assay would be used for further analysis
  tcell[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath,
    annotation = annotation
  )

  ### DNA accessibility data processing 
  DefaultAssay(tcell) <- "peaks"
  tcell <- FindTopFeatures(tcell, min.cutoff = 5)
  tcell <- RunTFIDF(tcell)
  tcell <- RunSVD(tcell)

  ############################################
  ### Data preprocessing of scRNA-seq data ###
  ############################################

  ### Set default assay to RNA
  DefaultAssay(tcell) <- "RNA"
  
  ### We normalized the gene expression through sctransform
  tcell <- SCTransform(tcell)  

  ### Dimensional reduction based on scRNA-seq data only
  tcell <- RunPCA(tcell)

  ### Lastly, we built a joint UMA visualization
  tcell <- RunUMAP(
    object = tcell,
    nn.name = "weighted.nn",
    assay = "RNA",
    verbose = TRUE
  )

  return(tcell)
}

### We processed the single-cell multiome data for each sample separately
naive_t <- Multiome_preprocessing("10XSC011-01/", "Naive_T")
active_t <- Multiome_preprocessing("10XSC011-02/", "Active_T")
exhaust_t72 <- Multiome_preprocessing("10XSC009-02/", "Exhausted_T_72hr")
exhaust_t96 <- Multiome_preprocessing("10XSC009-01/", "Exhausted_T_96hr")

### save each seurat object
saveRDS(naive_t, paste0(work_path, "naive_t.rds"))
saveRDS(active_t, paste0(work_path, "active_t.rds"))
saveRDS(exhaust_t72, paste0(work_path, "exhausted_t72.rds"))
saveRDS(exhaust_t96, paste0(work_path, "exhausted_t96.rds"))
