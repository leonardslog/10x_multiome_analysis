library(Signac)
library(Seurat)
library(rtracklayer)
library(GenomeInfoDb)
library(future)
plan(multisession, workers = 2) # parallelization
options(future.globals.maxSize= 11000*1024^2) # e.g. 6000gb*1024^2

workingdir <- ""
setwd(workingdir)

# load metadata 
meta <- read.csv(
  file = '2021-08-25_Marm038_brainstem_sl9_rxn1.atac.per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)

# create custom annotation
annotations <- rtracklayer::import('Callithrix_jacchus_cj1700_1/ncbi_dataset/data/GCF_009663435.1/genomic.gtf')
colnames(mcols(annotations))[colnames(mcols(annotations)) == "gene"] <- "gene_name" # gene_name col required for CreateChromatinAssay() annotation argument
annotations$tx_id <- annotations$transcript_id # tx_id possibly required for CoveraePlot downstream
annotations$gene_name <- annotations$gene # gene_name col required for CreateChromatinAssay() annotation argument
colnames(mcols(annotations))


# chromosome labels do not match so need to change them manually by vector mapping with UCSC Alias map:
alias_map_df <- read.csv("Callithrix_jacchus_cj1700_1/GCF_009663435.1.chromAlias.txt", 
                         header = TRUE, sep = "\t"
)
head(alias_map_df)
colnames(alias_map_df)[1] <- "refseq"
head(alias_map_df)

alias_map_vector <- setNames(alias_map_df$ucsc, alias_map_df$refseq) # try ucsc for downstream compatability
alias_map_vector <- alias_map_vector[!is.na(alias_map_vector) & alias_map_vector != ""] # remove NAs and blanks
head(alias_map_vector)

# Get current sequence names to change
current_seqlevels <- seqlevels(annotations)

# Map current names to new names using the vector (This will replace existing names if a match is found)
new_seqlevels <- alias_map_vector[as.character(current_seqlevels)]

# Assign new sequence names to the GRanges object (Note: this might drop ranges with NA seqnames if you have any)
seqlevels(annotations) <- new_seqlevels
annotations <- annotations[!is.na(annotations$gene_biotype)] # NAs will introduce error during TSSEnrichment
seqlevelsStyle(annotations) # now UCSC

######################################
# create multiome object and ATAC QC #
######################################

# load raw feaatures and remove droplets designated cell non-calls by CellRanger
counts <- Read10X_h5(filename = "2021-08-25_Marm038_brainstem_sl9_rxn1.atac.raw_feature_bc_matrix.h5")

multiome_data <- CreateSeuratObject(counts = counts$`Gene Expression`, meta.data = meta, assay = "RNA")
multiome_data <- subset(multiome_data, subset = is_cell == 1)

atac <- CreateSeuratObject(counts = counts$Peaks, meta.data = meta, assay = "ATAC")
atac <- subset(atac, subset = is_cell == 1)

multiome_data[["ATAC"]] <- CreateChromatinAssay(
  counts = LayerData(atac), 
  sep = c(":", "-"),
  genome = "calJac4",
  annotation = annotations,
  fragments = './2021-08-25_Marm038_brainstem_sl9_rxn1.atac_fragments.tsv.gz',
)
rm(atac, counts, alias_map_df, alias_map_vector, current_seqlevels, new_seqlevels)
gc()

# ATAC QC
DefaultAssay(multiome_data) <- "ATAC"
multiome_data <- NucleosomeSignal(object = multiome_data)
gc()
multiome_data <- TSSEnrichment(object = multiome_data, fast = FALSE)
gc()

# visualize above QC
DensityScatter(multiome_data, x="nCount_ATAC", y="TSS.enrichment", 
               log_x=TRUE, quantiles=TRUE)
VlnPlot(
  object = multiome_data,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = FALSE,
  ncol = 3,
)

# remove outliers
multiome_data_filtered <- subset(
  x = multiome_data,
  subset = nCount_ATAC < 74000 &
    nCount_ATAC > 2000 &
    TSS.enrichment > 2.7 &
    nucleosome_signal < 1.75
)
VlnPlot(
  object = multiome_data_filtered,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = FALSE,
  ncol = 3,
)
rm(multiome_data)
gc()


# normalize the gene expression and chromatin accessibility data
DefaultAssay(multiome_data_filtered) <- "RNA"
multiome_data_filtered <- NormalizeData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- FindVariableFeatures(multiome_data_filtered, assay="RNA", nfeatures = 5000)
multiome_data_filtered <- ScaleData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- RunPCA(multiome_data_filtered)

DefaultAssay(multiome_data_filtered) <- "ATAC"
multiome_data_filtered <- FindTopFeatures(multiome_data_filtered, min.cutoff = 5)
multiome_data_filtered <- RunTFIDF(multiome_data_filtered)
multiome_data_filtered <- RunSVD(multiome_data_filtered)

# checkpoint save
SaveSeuratRds(object = multiome_data_filtered, file = "multiome_data_filtered.Rds")

# save .h5ad for MapMyCells
library(scCustomize)
library(reticulate)

use_condaenv("single-cell") 
as.anndata(
  multiome_data_filtered,
  file_path = "",
  file_name= "multiome_data_filtered.h5ad",
  assay = "RNA",
  main_layer = "counts",
  other_layers = NULL,
  transer_dimreduc = FALSE,
  verbose = TRUE,
)