library(Seurat)
library(presto)
library(Signac)
library(GenomicRanges)
library(future)
plan(multisession, workers = 2) # parallelization
options(future.globals.maxSize= 11000*1024^2) # e.g. 6000gb*1024^2

workingdir <- ""
setwd(workingdir)

# load mapmycell annotations
mapmycells_WMB <- read.csv(file = "multiome_data_filtered_10xWholeMouseBrain(CCN20230722)/multiome_data_filtered_10xWholeMouseBrain(CCN20230722).csv",
                             skip = 4, header = TRUE)
head(mapmycells_WMB)
# set cell_id as rownames
rownames(mapmycells_WMB) <- mapmycells_WMB$cell_id

data.frame(Cell_counts=head(sort(table(mapmycells_WMB$class_name),decreasing=T),34))
# most well represented cells
# Cell_counts.Var1 Cell_counts.Freq
# 1       31 OPC-Oligo            10258
# 2          23 P Glut             2550
# 3      30 Astro-Epen             1407
# 4          34 Immune              861
# 5         20 MB GABA              269
# 6         24 MY Glut              222
# 7         19 MB Glut              155
# 8          26 P GABA              134
# 9        33 Vascular               67
# 10     01 IT-ET Glut               62
# 11        21 MB Dopa               60
# 12    05 OB-IMN GABA               54
# 13   13 CNU-HYa Glut               22
# 14        12 HY GABA               20
# 15        14 HY Glut               14
# 16        29 CB Glut               13
# 17    04 DG-IMN Glut                9
# 18   11 CNU-HYa GABA                9
# 19 02 NP-CT-L6b Glut                8
# 20   06 CTX-CGE GABA                8
# 21     22 MB-HB Sero                7
# 22        27 MY GABA                5
# 23   07 CTX-MGE GABA                4
# 24     17 MH-LH Glut                3
# 25        28 CB GABA                2
# 26            32 OEC                2
# 27     03 OB-CR Glut                1
# 28   08 CNU-MGE GABA                1

# load multiome object and annotate w/ mapmycells
data <- LoadSeuratRds("01_multiome_data_filtered.Rds")

# sanity checks:
# check cell_id/barcode equivalence btwn datasets and btwn mapmycells and data
head(rownames(mapmycells_WMB)) == head(colnames(data[["ATAC"]])) # TRUE
head(rownames(mapmycells_WMB)) == head(colnames(data[["RNA"]])) # TRUE
df <- (colnames(data[["RNA"]]) == colnames(data[["ATAC"]]))
FALSE %in% df # RNA and ATAC columns are ordered the same
rm(df)

## if data and annotations aren't in same order, comment out following:
## index to reorder new_annotations (match(A, B) returns positions of A elements in B)
# match_index <- match(colnames(data[["RNA"]]), rownames(mapmycells_WMB))

## reorder the cell annotations data frame
# mapmycells_mouse_reordered <- mapmycells_mouse[match_index, ]

# all(rownames(mapmycells_mouse_reordered) == colnames(rna)) # sanity check: TRUE
# all(rownames(mapmycells_mouse_reordered) == colnames(atac)) # sanity check: TRUE
# rm(mapmycells_mouse)

# add MapMyCells metadata
data <- AddMetaData(object = data, metadata = mapmycells_WMB)

# inspect dataset
DefaultAssay(data) <- "RNA"
print(data[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  CACNA1C, CCSER1, LOC100393071, DLGAP1, CNTNAP5 
## Negative:  LOC118151440, CHN2, SEMA3E, IL1RAP, LOC100385822 
## PC_ 2 
## Positive:  GLI3, GLIS3, RFX4, ADCY2, ITPRID1 
## Negative:  CDH8, HS3ST5, SGK1, OSBPL6, THSD7B 
## PC_ 3 
## Positive:  LOC103788313, ARHGAP15, MEF2C, DOCK8, INPP5D 
## Negative:  ADGRB3, GRID2, ERBB4, NEGR1, ANKS1B 
## PC_ 4 
## Positive:  PDGFRA, LHFPL3, LOC118147387, MMP16, TMEM132D 
## Negative:  SLC4A4, MGAT4C, GASK1A, BMPR1B, GLI3 
## PC_ 5 
## Positive:  FOXP2, PTCHD4, GAD2, COL25A1, LOC118155279 
## Negative:  ANO3, ARPP21, DGKG, TNR, NFIB 

VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca") + NoLegend()
par(mar = c(5, 4, 4, 2) + 0.1)
ElbowPlot(data) # using dims 1-11

# link assays / find clusters
data <- FindNeighbors(data, dims = 1:11)
data <- FindClusters(data, resolution = 1.5)
data <- RunUMAP(data, dims = 1:11)
DimPlot(data, reduction = "umap")
DimPlot(data, reduction = "umap", group.by = "class_name", repel = T, label = T)

# DE analysis-------------------------------------------------------------------
library(ggplot2)

# find DE features between "23 P Glut" and "30 Astro-Epen"
DefaultAssay(data) <- "RNA"
Idents(data) <- "class_name"

oligo_pglut.de.markers <- FindMarkers(data, ident.1 = "31 OPC-Oligo", ident.2 = "23 P Glut")
head(oligo_pglut.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# EPHA7            0  -4.990820 0.043 0.962         0
# ADCY8            0  -5.121001 0.052 0.971         0
# SRRM3            0  -4.881668 0.032 0.945         0
# LOC103787440     0  -5.257025 0.028 0.939         0
# TENM4            0  -4.886707 0.085 0.988         0
# LOC103796371     0  -4.610896 0.056 0.959         0

range(oligo_pglut.de.markers$avg_log2FC)
# [1] -8.254263  5.958942
write.csv(oligo_pglut.de.markers,file = "oligo_pglut_de_markers.csv")


# DA analysis-------------------------------------------------------------------
library(Biostrings)
genome <- readDNAStringSet("Callithrix_jacchus_cj1700_1/ncbi_dataset/data/GCF_009663435.1/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.fna")
genome.names <- names(genome)
names(genome)[1:10]

alias_map_df <- read.csv("Callithrix_jacchus_cj1700_1/GCF_009663435.1.chromAlias.txt", 
                         header = TRUE, sep = "\t"
)

# name formatting: split the names of the genome by " " and take the first index
name_index <- sub(" .*", "", genome.names)
names(genome) <- name_index

head(alias_map_df)
colnames(alias_map_df)[1] <- "refseq"
head(alias_map_df)
alias_map_vector <- setNames(alias_map_df$ucsc, alias_map_df$refseq)

# Get current sequence names to change
current_seqnames <- names(genome)

# Map current names to new names using the vector
new_seqnames <- alias_map_vector[as.character(current_seqnames)]
names(genome) <- new_seqnames
GenomeInfoDb::seqlevelsStyle(genome) # UCSC

# this will take a minute without 'presto' package
gc()
DefaultAssay(data) <- "ATAC"
da_peaks_opc_pglut <- FindMarkers(
  object = data,
  ident.1 = "31 OPC-Oligo",
  ident.2 = "23 P Glut",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks_opc_pglut)
## p_val avg_log2FC pct.1 pct.2 p_val_adj
## chr1-84300555-84301446       0  -5.051686 0.048 0.685         0
## chr6-150538031-150538933     0   3.210641 0.615 0.096         0
## chr4-39431736-39432775       0  -5.919404 0.013 0.513         0
## chr2-201744178-201745109     0  -5.201473 0.020 0.507         0
## chr5-44216315-44217194       0  -3.370630 0.079 0.562         0
## chr7-65237582-65238407       0  -6.923535 0.006 0.488         0

range(da_peaks_opc_pglut$avg_log2FC)
# [1] -11.252793   3.936752

plot1 <- VlnPlot(
  object = data,
  features = rownames(da_peaks_opc_pglut)[1],
  pt.size = 0.1,
  idents = c("31 OPC-Oligo","23 P Glut")
)
plot2 <- FeaturePlot(
  object = data,
  features = rownames(da_peaks_opc_pglut)[1],
  pt.size = 0.1
)

plot1 | plot2

# subset open regions
open_oligo <- rownames(da_peaks_opc_pglut[da_peaks_opc_pglut$avg_log2FC > 2.5, ])
open_pglut <- rownames(da_peaks_opc_pglut[da_peaks_opc_pglut$avg_log2FC < -5, ])

closest_genes_oligo <- ClosestFeature(data, regions = open_oligo)
# Warning message:
# In ClosestFeature(data, regions = open_oligo) :
# The following seqlevels present in query regions are not present
# in the supplied gene annotations and will be removed: chrUn.000463F.qpd.obj, chrUn.000610F.qpd.obj, 
# chrUn.001016F.qpd.obj, Super.Scaffold.100049, Super.Scaffold.100088

closest_genes_pglut <- ClosestFeature(data, regions = open_pglut)
# Warning message:
# In ClosestFeature(data, regions = open_pglut) :
# The following seqlevels present in query regions are not present
# in the supplied gene annotations and will be removed: chrUn.000149F.qpd.subseq.3300941, 
# chrUn.000601F.qpd.obj, chrUn.000722F.qpd.obj, chrUn.000738F.qpd.obj, Super.Scaffold.100060, 
# Super.Scaffold.100071, Super.Scaffold.100092

head(closest_genes_oligo)
# source type score phase gene_id transcript_id          db_xref gbkey   gene
# 1 Gnomon gene    NA    NA  B3GNT7               GeneID:100394256  Gene B3GNT7
# 2 Gnomon gene    NA    NA   ARNTL               GeneID:100396695  Gene  ARNTL
# 3 Gnomon gene    NA    NA   PSMC4               GeneID:100387347  Gene  PSMC4
# 4 Gnomon gene    NA    NA  NKX2-2               GeneID:100389271  Gene NKX2-2
# 5 Gnomon gene    NA    NA    FA2H               GeneID:100387942  Gene   FA2H
# 6 Gnomon gene    NA    NA    HIP1               GeneID:100404703  Gene   HIP1

head(closest_genes_pglut)
# source type score phase      gene_id transcript_id          db_xref gbkey         gene
# 1 Gnomon gene    NA    NA        RPS10               GeneID:103792402  Gene        RPS10
# 2 Gnomon gene    NA    NA LOC118143226               GeneID:118143226  Gene LOC118143226
# 3 Gnomon gene    NA    NA LOC100406071               GeneID:100406071  Gene LOC100406071
# 4 Gnomon gene    NA    NA LOC118148948               GeneID:118148948  Gene LOC118148948
# 5 Gnomon gene    NA    NA        NDRG4               GeneID:100389623  Gene        NDRG4
# 6 Gnomon gene    NA    NA       MTMR14               GeneID:100397575  Gene       MTMR14

# Region stats (GC content, etc) 
data <- RegionStats(object = data, genome = genome)

# find DA peaks overlapping gene of interest

# limit plot to specific celltypes
idents.plot <- c("31 OPC-Oligo", "23 P Glut", "30 Astro-Epen", "34 Immune")

data <- SortIdents(data)

regions_highlight <- subsetByOverlaps(StringToGRanges(open_oligo), LookupGeneCoords(data, "B3GNT7"))
regions_highlight

# redefined annotation metadata to accommodate inconsistency in requirements from TSSEnrichment to CoveragePlot
annotations <- rtracklayer::import('Callithrix_jacchus_cj1700_1/ncbi_dataset/data/GCF_009663435.1/genomic.gtf')
current_seqlevels <- seqlevels(annotations)
new_seqlevels <- alias_map_vector[as.character(current_seqlevels)]
seqlevels(annotations) <- new_seqlevels
data[["ATAC"]]@annotation$gene_name <- data[["ATAC"]]@annotation$gene
data[["ATAC"]]@annotation$tx_id <- data[["ATAC"]]@annotation$transcript_id 
seqlevelsStyle(genome)
seqlevelsStyle(annotations)

# visualize genome tracks for selected genes
CoveragePlot(
  object = data,
  assay = "ATAC",
  expression.assay = "RNA",
  region = "B3GNT7",
  region.highlight = regions_highlight,
  extend.upstream = 2000,
  extend.downstream = 5000,
  idents = idents.plot,
  annotation = TRUE,
)

# compute gene activities
gene.activities <- GeneActivity(data)
# Warning message:
# In SingleFeatureMatrix(fragment = fragments[[x]], features = features,  :
# 219 features are on seqnames not present in the fragment file. These will be removed.

data[["GENE_ACTIVITIES"]] <- CreateAssayObject(counts = gene.activities)
data <- NormalizeData(
  object = data,
  assay = 'GENE_ACTIVITIES',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

SaveSeuratRds(object = data, file = "02_multiome_data_filtered.Rds")
