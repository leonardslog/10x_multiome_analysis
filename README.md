# Processing and Analysis of 10X Multiomic data

## Overview

To familiarize myself with RNA and ATAC-seq analysis workflows, I opted to work raw data publicly available via the Allen Institute's [Brain Knowledge Platform](https://brain-map.org/bkp). These analyses were to be carried out on a laptop with 8 cores and 64GB of memory, so to address computational resource limitations, I analyzed 10X Multiomic data generated from brainstem tissue in the common marmoset *Callithrix jacchus*.

Much of the following workflow and code are adapted from the following tutorials/vignettes: \
https://stuartlab.org/signac/articles/overview \
https://satijalab.org/seurat/articles/get_started_v5_new \
https://ngs101.com/tutorials/ \
https://github.com/mousepixels/sanbomics_scripts

## System setup

`Seurat` and `Signac` require `future` for parallelization. I'm performing these analyses on a laptop with 8 cores and 64GB, and during some steps had to adjust my memory allocations to accomodate the processes. `rtracklayer` and `GenomeInfoDb` are needed for constructing the annotation object. At the time of this exercise I used the developer version of `Seurat` to grab the fix for a [plotting bug](https://github.com/satijalab/seurat/issues/10180) I encountered in downstream plotting. `scCustomize` and `reticulate` are required to save the Seurat object to be constructed as an .h5ad object required for MapMyCells. `reticulate` requires a separate python installation that has the `anndata` module installed. 

~~~
library(Seurat) # v5.3.1.9001
library(Signac)
library(rtracklayer)
library(GenomeInfoDb)
library(future)
plan(multisession, workers = 2) # parallelization
options(future.globals.maxSize= 11000*1024^2) # e.g. 6000gb*1024^2
~~~

## Pre-processing

Loading metadata for the raw counts:
```
meta <- read.csv(
  file = '2021-08-25_Marm038_brainstem_sl9_rxn1.atac.per_barcode_metrics.csv',
  header = TRUE,
  row.names = 1)
```

For the custom annotation, I have to use the same _C. jacchus_ genome build as what the reads were aligned to (NCBI RefSeq Assembly ID: [GCF_009663435.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009663435.1/), or 'caljac4' in the UCSC browser.
```
annotations <- rtracklayer::import('Callithrix_jacchus_cj1700_1/ncbi_dataset/data/GCF_009663435.1/genomic.gtf')
```
Downstream, I learned some functions expect certain column names in the Assay objects I'll create shortly and will call errors without them, so they're added here.
```
colnames(mcols(annotations))[colnames(mcols(annotations)) == "gene"] <- "gene_name" # gene_name required for CreateChromatinAssay() annotation parameter
annotations$tx_id <- annotations$transcript_id # tx_id column potentially required to prevent CoveragePlot() error
annotations$gene_name <- annotations$gene # gene_name col required for CreateChromatinAssay() annotation argument
colnames(mcols(annotations))
```

Seurat's CreateChromatinAssay() also requires UCSC format for the chromosome labels of the reference (currently RefSeq), so I changed them by vector mapping with the Alias.txt downloaded from UCSC.
```
alias_map_df <- read.csv("Callithrix_jacchus_cj1700_1/GCF_009663435.1.chromAlias.txt", 
                         header = TRUE, sep = "\t"
)
head(alias_map_df)
colnames(alias_map_df)[1] <- "refseq"
head(alias_map_df)
       refseq assembly    genbank ncbi ucsc
1 NC_025586.1       MT              MT chrM
2 NC_048383.1     chr1 CM018917.1    1 chr1
3 NC_048384.1     chr2 CM018918.1    2 chr2
4 NC_048385.1     chr3 CM018919.1    3 chr3
5 NC_048386.1     chr4 CM018920.1    4 chr4
6 NC_048387.1     chr5 CM018921.1    5 chr5

alias_map_vector <- setNames(alias_map_df$ucsc, alias_map_df$refseq)
alias_map_vector <- alias_map_vector[!is.na(alias_map_vector) & alias_map_vector != ""] # remove NAs and blanks in case present
head(alias_map_vector)
NC_025586.1 NC_048383.1 NC_048384.1 NC_048385.1 NC_048386.1 NC_048387.1 
     "chrM"      "chr1"      "chr2"      "chr3"      "chr4"      "chr5" 
```

I then constructed a vector of the current chromosomes/scaffolds in the genome that need to be changed, mapped them to the names to names to UCSC names, and assigned them to the genome/annotation object.
```
current_seqlevels <- seqlevels(annotations)
new_seqlevels <- alias_map_vector[as.character(current_seqlevels)]
seqlevels(annotations) <- new_seqlevels
```

Another retroactive fix to a downstream error involves the presence of NA's in the gene_biotype column of the annotation.
```
annotations <- annotations[!is.na(annotations$gene_biotype)] # NAs will introduce error during TSSEnrichment
seqlevelsStyle(annotations) # confirmed that new style is UCSC
```

## Create Multiome Seurat object

Counts were available as a raw feature matrix, which after loading I filtered to remove all droplets deemed non-calls by the CellRanger pipeline.

```
counts <- Read10X_h5(filename = "2021-08-25_Marm038_brainstem_sl9_rxn1.atac.raw_feature_bc_matrix.h5")

multiome_data <- CreateSeuratObject(counts = counts$`Gene Expression`, meta.data = meta, assay = "RNA")
multiome_data <- subset(multiome_data, subset = is_cell == 1)

atac <- CreateSeuratObject(counts = counts$Peaks, meta.data = meta, assay = "ATAC")
atac <- subset(atac, subset = is_cell == 1)
```

Next. I created the chromatin ASSAY with the ATAC feature matrix, which needs both the fragments and the annotations object constructed above. I also removed the objects no longer needed to free up memory.
```
multiome_data[["ATAC"]] <- CreateChromatinAssay(
  counts = LayerData(atac), 
  sep = c(":", "-"),
  genome = "calJac4",
  annotation = annotations,
  fragments = './2021-08-25_Marm038_brainstem_sl9_rxn1.atac_fragments.tsv.gz',
)
rm(atac, counts, current_seqlevels, new_seqlevels)
gc()
```
Per-cell uality control metrics were calculated and visualized to inform additional filtering.
```
DefaultAssay(multiome_data) <- "ATAC"
multiome_data <- NucleosomeSignal(object = multiome_data)
multiome_data <- TSSEnrichment(object = multiome_data, fast = FALSE)

DensityScatter(multiome_data, x="nCount_ATAC", y="TSS.enrichment", 
               log_x=TRUE, quantiles=TRUE)
VlnPlot(
  object = multiome_data,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  pt.size = FALSE,
  ncol = 3
)
```

I filtered cells with cuttoffs for counts, TSS enrichment, and nucleosome signal based on the bottom and top quantiles shown in the density scatterplot.

```
multiome_data_filtered <- subset(
  x = multiome_data,
  subset = nCount_ATAC < 74000 &
    nCount_ATAC > 2000 &
    TSS.enrichment > 2.7 &
    nucleosome_signal < 1.75
)
rm(multiome_data)
gc()
```
## Normalization
**_Note_**: MapMyCells does not require normalization for cell annotation, and developers urge caution in interpreting results from log-normalized data, so if performing this step first, make sure to point to the right Seurat object layer when formatting input for it.

```
DefaultAssay(multiome_data_filtered) <- "RNA"
multiome_data_filtered <- NormalizeData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- FindVariableFeatures(multiome_data_filtered, assay="RNA", nfeatures = 5000)
multiome_data_filtered <- ScaleData(multiome_data_filtered, assay="RNA")
multiome_data_filtered <- RunPCA(multiome_data_filtered)

DefaultAssay(multiome_data_filtered) <- "ATAC"
multiome_data_filtered <- FindTopFeatures(multiome_data_filtered, min.cutoff = 5)
multiome_data_filtered <- RunTFIDF(multiome_data_filtered)
multiome_data_filtered <- RunSVD(multiome_data_filtered)
```
MapMyCells requires an anndata (.h5ad) or .csv for input, so converting (and saving an R object for post-annotation analyses).

```
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

SaveSeuratRds(object = multiome_data_filtered, file = "multiome_data_filtered.Rds")
```

## Cell Annotation

After annotating the cells using the MapMyCells browser platform using the 10x Whole Mouse Brain Taxonomy (CCN20230722) as the reference. I restarted a new environment and added the output as metadata to the previously saved Seurat object. Rownames renamed as the cell_ID to match the metadata of the Seurat object.

```
library(Seurat)
library(presto)
library(Signac)
library(GenomicRanges)
library(future)
plan(multisession, workers = 2) # parallelization
options(future.globals.maxSize= 11000*1024^2) # e.g. 6000gb*1024^2

mapmycells_WMB <- read.csv(file = "multiome_data_filtered_10xWholeMouseBrain(CCN20230722)/multiome_data_filtered_10xWholeMouseBrain(CCN20230722).csv",
                             skip = 4, header = TRUE)
colnames(mapmycells_WMB)
# [1] "cell_id"                             "class_label"                        
# [3] "class_name"                          "class_bootstrapping_probability"    
# [5] "subclass_label"                      "subclass_name"                      
# [7] "subclass_bootstrapping_probability"  "supertype_label"                    
# [9] "supertype_name"                      "supertype_bootstrapping_probability"
#[11] "cluster_label"                       "cluster_name"                       
#[13] "cluster_alias"                       "cluster_bootstrapping_probability" 

rownames(mapmycells_WMB) <- mapmycells_WMB$cell_id
```
From the distibution of celltype counts (N=16227), over 60% of the cells are identified as oligodendrocytes or oligodendrocyte precursor cells, with the 2nd and 3rd largest clusters consisting of Pons and glutamatergic neuron cells (#23) and astrocytes or ependymal cells (#30).
```
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
```



## Differential Expression analysis


## Differential Accessibility analysis


## Celltype prediction with primate data (_Microcebus murinus_) and SingleR
I downloaded RNA-seq data derived from brainstem and cortex tissue from the [Tabula murnius]("https://tabula-microcebus.sf.czbiohub.org/whereisthedata") and converted the .h5ad to a Seurat object.
```
library(anndataR)
library(Seurat)
library(SingleR)

# reference data preparation
reference <- read_h5ad("~/Downloads/microcebus_murinus_brainstem_and_cortex.h5ad")
reference <- reference$as_Seurat()
```

Similar to the annotation construction and chromosome renaming prior to generating the Multiome Seurat object, here I needed to convert the reference's rownames from ENSEMBL gene identifiers to gene names to match the query data. This was accomplished by querying the ENSEMBL database with ```biomaRt```.
```
library(biomaRt)
listMarts() 
#               biomart                version
#1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 115
#2   ENSEMBL_MART_MOUSE      Mouse strains 115
#3     ENSEMBL_MART_SNP  Ensembl Variation 115
#4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 115

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
mouselemur_mart <- useDataset(dataset = 'mmurinus_gene_ensembl', mart = ensembl, verbose = TRUE)
attributes <- c("ensembl_gene_id", "external_gene_name")
gene_id_data <- getBM(attributes = attributes, mart = mouselemur_mart)
head(gene_id_data)

id_to_symbol <- setNames(gene_id_data$external_gene_name, gene_id_data$ensembl_gene_id)

# get ref's counts matrix and assign new rownames
rownames(reference) <- id_to_symbol[rownames(reference)]
reference <- subset(reference, features = rownames(reference)[rownames(reference) != ""])
```
I ran SingleR under the Wilcoxon Rank Sum Test, first using log-normalized data. The result is a large dataframe, which was saved as an R object.
```
data <- LoadSeuratRds("~/Documents/shixuan_liu/data/marmoset_brainstem/multiome_data_final.Rds")
norm_counts <- LayerData(data, assay = "RNA", layer = 'data') # LogNormalized
raw_counts <- LayerData(data, assay = "RNA", layer = 'counts')

# run SingleR
cellannotation_microcebus <- SingleR(test = norm_counts,
                  ref = reference[["RNA"]]$X, 
                  labels = reference$cell_type,
                  de.method = 'wilcox')
saveRDS(cellannotation_microcebus,file = "singleR_microcebus_cellannotation.Rds")
```
I visualized cell x gene expression patterns were visualized with ```pheatmap```.
```
library(pheatmap)
plotScoreHeatmap(cellannotation_microcebus)
plotDeltaDistribution(cellannotation_microcebus, ncol = 4, dots.on.top = FALSE)
```
Finally, I added these annotations to the original (and now very large) Seurat object's metadata to faciliate comparison of celltype annotations.
```
rownames(cellannotation_microcebus)[1:5] # make sure you have cell IDs
data <- Seurat::AddMetaData(data, cellannotation_microcebus$pruned.labels, col.name = 'SingleR_annt_microcebus_ref')
saveRDS(data,file = "~/Documents/shixuan_liu/data/marmoset_brainstem/multiome_data_final_v2.Rds")
```
As one might expect, the majority of the cells are identified as oligodendroctyes, but 38% of the cells are assigned to either neurons (not subdivided as in MapMyCells) or astrocytes. Thee MapMyCells WMB annotations and the SingleR Microcebus annotations do not have perfectly congruent celltypes, but between both strategies there is high similarity in number of assignments to either oligodendrocytes or OPCs (MMC-WMB = 9098, SingleR-Mm = 10258), as well as the labeling of either astrocytes or ependymal cells (MMC-WMB = 1407, SingleR-Mm = 1588). Furthermore, 303 cells were also pruned during the SingleR prediction due to low quality annotation.

```
microcebus_annotations <- readRDS("Documents/shixuan_liu/data/marmoset_brainstem/singleR_microcebus_cellannotation.Rds")
data <- readRDS("Documents/shixuan_liu/data/marmoset_brainstem/multiome_data_final_v2.Rds")
data <- SetIdent(data, value = 'SingleR_annt_microcebus_ref')
DimPlot(data, reduction = "umap", repel = T, label = T) #+ NoLegend()
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
```
Visualization of shared gene features between the _Microcebus murinus_ and _Callithrix jacchus_ references, and Multiome data show there are 860 genes in the _Microcebus_ genome that are absent in the _Callithrix_. While irrelevant for this experiment due their absence in a dataset derived from a subset of brain tissue with very specific expression profiles, suboptimal genome selection for annotation can drastically impact downstream results.
```
library(VennDiagram)
venn.diagram(
  disable.logging = TRUE,
  margin =.1,
  x = list(microcebus_genes, callithtrix_genes, rna_genes, gene_activities),
  category.names = c("Microcebus" , "Callithrix", "RNA Assay", "Gene Activities Assay"),
  filename = 'microcebus_vs_callithrix_venn_diagramm_multiassay.png',
  output=TRUE
)
```

## References

### Data

Dataset: Fenna M. Krienen, Kirsten M. Levandowski, Heather Zaniewski, Ricardo C.H. del Rosario, Margaret E. Schroeder, Melissa Goldman, Alyssa Lutservitz, Qiangge Zhang, Katelyn X. Li, Victoria F. Beja-Glasser, Jitendra Sharma, Tay Won Shin, Abigail Mauermann, Alec Wysoker, James Nemesh, Seva Kashin, Josselyn Vergara, Gabriele Chelini, Jordane Dimidschstein, Sabina Berretta, Ed Boyden, Steven A. McCarroll, Guoping Feng (2023) A Marmoset Brain Cell Census Reveals Persistent Influence of Developmental Origin on Neurons. Available from https://assets.nemoarchive.org/dat-1je0mn3

UCSC Genome Browser assembly ID: "calJac4". Perez et al. The UCSC Genome Browser database: 2025 update. Nucleic Acids Research 2025 PMID: 39460617, DOI: 10.1093/nar/gkae974 

### Literature

Ezran, Camille, et al. "Mouse lemur cell atlas informs primate genes, physiology and disease." Nature 644.8075 (2025): 185-196. 

Gontarz, Paul, et al. "Comparison of differential accessibility analysis strategies for ATAC-seq data." Scientific reports 10.1 (2020): 10150.

Krienen, Fenna M., et al. "A marmoset brain cell census reveals regional specialization of cellular identities." Science advances 9.41 (2023): eadk3986. 

Teo, Alan Yue Yang, et al. "Best practices for differential accessibility analysis in single-cell epigenomics." Nature communications 15.1 (2024): 8805. 

Yan, Feng, et al. "From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis." Genome biology 21.1 (2020): 22. 

Yao, Zizhen, et al. "A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain." Nature 624.7991 (2023): 317-332. 

