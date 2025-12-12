library(anndataR)
library(Seurat)
library(SingleR)


# reference data preparation
reference <- read_h5ad("microcebus_murinus_brainstem_and_cortex.h5ad")
ref <- reference$as_Seurat()
rm(reference)

# need to convert rownames from gene identifiers to gene names
library(biomaRt)
listMarts() 
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl) 
mouselemur_mart <- useDataset(dataset = 'mmurinus_gene_ensembl', mart = ensembl, verbose = TRUE)
attributes <- c("ensembl_gene_id", "external_gene_name")
gene_id_data <- getBM(attributes = attributes, mart = mouselemur_mart)
head(gene_id_data)

id_to_symbol <- setNames(gene_id_data$external_gene_name, gene_id_data$ensembl_gene_id)

# get ref's counts matrix and assign new rownames
ref_copy <- ref
rownames(ref_copy) <- id_to_symbol[rownames(ref_copy)]
ref_copy <- subset(ref_copy, features = rownames(ref_copy)[rownames(ref_copy) != ""])

# query data
data <- LoadSeuratRds("03_multiome_data_filtered.Rds")
norm_counts <- LayerData(data, assay = "RNA", layer = 'data') # LogNormalized
raw_counts <- LayerData(data, assay = "RNA", layer = 'counts')

# run SingleR
cellannotation_microcebus <- SingleR(test = norm_counts,
                  ref = ref_copy[["RNA"]]$X, 
                  labels = ref_copy$cell_type,
                  de.method = 'wilcox')
saveRDS(cellannotation_microcebus,file = "04_singleR_microcebus_cellannotation.Rds")

# Inspect quality of the predictions
library(pheatmap)
plotScoreHeatmap(cellannotation_microcebus)
plotDeltaDistribution(cellannotation_microcebus, ncol = 4, dots.on.top = FALSE)

# add to final seurat to save
rownames(cellannotation_microcebus)[1:5] # make sure you have cell IDs
data <- Seurat::AddMetaData(data, cellannotation_microcebus$pruned.labels, col.name = 'SingleR_annt_microcebus_ref')
saveRDS(data,file = "04_multiome_data_final.Rds")

# Visualize UMAP
library(Seurat)
microcebus_annotations <- readRDS("singleR_microcebus_cellannotation.Rds")
data <- readRDS("04_multiome_data_final.Rds")
data <- SetIdent(data, value = 'SingleR_annt_microcebus_ref')

DimPlot(data, reduction = "umap", repel = T, label = T)
# Error in data.frame(..., check.names = FALSE) : 
#   arguments imply differing number of rows: 16227, 15924
# In addition: Warning message:
#   Removing 303 cells missing data for vars requested 
# referenced here: https://github.com/satijalab/seurat/issues/10180

data.frame(Cell_counts=head(sort(table(microcebus_annotatoins$labels),decreasing=T),34))
# Cell_counts.Var1 Cell_counts.Freq
# 1                 oligodendrocyte             8202
# 2                          neuron             4533
# 3                       astrocyte             1578
# 4  oligodendrocyte precursor cell              896
# 5                GABAergic neuron              351
# 6                      macrophage              249
# 7            glutamatergic neuron              183
# 8                    stromal cell               91
# 9      capillary endothelial cell               30
# 10       myelinating Schwann cell               21
# 11   non-myelinating Schwann cell               18
# 12                         B cell               14
# 13                     lymphocyte               13
# 14                       pericyte               13
# 15                 ependymal cell               10
# 16 choroid plexus epithelial cell                9
# 17            leptomeningeal cell                8
# 18                        unknown                8