library(Seurat)
library(presto)
library(Signac)
# library(GenomicRanges)

workingdir <- ""
setwd(workingdir)

# reload data
data <- LoadSeuratRds("02_multiome_data_filtered.Rds")

# CZBiohub marker genes for 4 most represented MapMyCells celltypes
all_markers <- c("ST18", "CTNNA3", "RNF220", "PIP4K2A", "MBP", "TMEM144", "PDE4B", 
             "PHLPP1", "EDIL3", "PRUNE2", "ELMO1", "TTLL7", "C10orf90", "ENPP2", 
             "DOCK5", "MAP7", "MAN2A1", "QKI", "TMTC2", "SHTN1", "PLCL1", "PLEKHH1", 
             "TF", "PXK", "SLC24A2", "PTPRZ1", "TNR", "LHFPL3", "PCDH15", "EPN2", 
             "SOX6", "SCD5", "CA10", "LUZP2", "VCAN", "LRRC4C", "NXPH1", "NOVA1", "XYLT1", 
             "MMP16", "KAT2B", "MEGF11", "SEMA5A", "SMOC1", "SLC35F1", "GALNT13", "AGAP1", 
             "DGKG", "CDH20", "ARPP21", "DLGAP2", "RALYL", "R3HDM1", "KHDRBS2", "NELL2", "SH3GL2", 
             "SATB2", "KHDRBS3", "SV2B", "HECW1", "GRIN2A", "KCNQ5", "RYR2", "CHRM3", "KCNIP4", "LDB2", 
             "PTPRD", "PRKCB", "IQCJ-SCHIP1", "CHN1", "PHACTR1", "CNKSR2", "TAFA1", "FSTL4",
             "SLC1A2", "NPAS3", "DTNA", "GPC5", "PITPNC1", "GPM6A", "BMPR1B", "RFX4", "MSI2",
             "ADGRV1", "GLIS3", "RYR3", "SLC4A4", "TRPS1", "ATP13A4", "ATP1A2", "ADCY2", "NTM", "ZNRF3", 
             "AQP4", "OBI1-AS1", "LRIG1", "PREX2")
oligo_markers <- c("ST18", "CTNNA3", "RNF220", "PIP4K2A", "MBP", "TMEM144", "PDE4B", 
                   "PHLPP1", "EDIL3", "PRUNE2", "ELMO1", "TTLL7", "C10orf90", "ENPP2", 
                   "DOCK5", "MAP7", "MAN2A1", "QKI", "TMTC2", "SHTN1", "PLCL1", "PLEKHH1", 
                   "TF", "PXK", "SLC24A2")
opc_markers <- c("PTPRZ1", "TNR", "LHFPL3", "PCDH15", "EPN2", "SOX6", "SCD5", 
                 "CA10", "LUZP2", "VCAN", "LRRC4C", "NXPH1", "NOVA1", "XYLT1", 
                 "MMP16", "PDE4B", "KAT2B", "MEGF11", "SEMA5A", "SMOC1", "SLC35F1", 
                 "GALNT13", "AGAP1", "DGKG", "CDH20")

pglut_markers <- c("ARPP21", "DLGAP2", "RALYL", "R3HDM1", "KHDRBS2", "NELL2", 
                   "SH3GL2", "SATB2", "KHDRBS3", "SV2B", "HECW1", "GRIN2A", "KCNQ5", 
                   "RYR2", "CHRM3", "KCNIP4", "LDB2", "PTPRD", "PRKCB", "IQCJ-SCHIP1", 
                   "CHN1", "PHACTR1", "CNKSR2", "TAFA1", "FSTL4")
astro_markers <- c("CDH20", "SLC1A2", "NPAS3", "DTNA", "GPC5", "PITPNC1", "PTPRZ1", 
                   "GPM6A", "BMPR1B", "RFX4", "MSI2", "ADGRV1", "GLIS3", "RYR3", 
                   "SLC4A4", "TRPS1", "ATP13A4", "ATP1A2", "ADCY2", "NTM", "ZNRF3", 
                   "AQP4", "OBI1-AS1", "LRIG1", "PREX2")

library(VennDiagram)
venn.diagram(
  disable.logging = TRUE,
  margin =.1,
  x = list(oligo_markers, opc_markers, pglut_markers, astro_markers),
  category.names = c("Oligodendrocytes" , "OPCs" , "Glutamatergic neurons", "Astrocytes"),
  filename = 'venn_diagramm.png',
  output=TRUE
)
intersect(oligo_markers, opc_markers) # "PDE4B"
intersect(astro_markers, opc_markers) # "CDH20"  "PTPRZ1"
all_markers <-all_markers[!all_markers %in% c("PDE4B", "CDH20",  "PTPRZ1")]

# link peaks to gene activity in all cell markers-------------------------------

data <- LinkPeaks(
  object = data,
  peak.assay = "ATAC",
  expression.assay = "GENE_ACTIVITIES",
  peak.slot = "counts",
  method = "pearson",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  genes.use = all_markers,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)

# differential expression-------------------------------------------------------
library(ggplot2)

DefaultAssay(data) <- "RNA"
Idents(data) <- "class_name"

# oligodendrocytes vs glutamatergic neurons
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
write.csv(oligo_pglut.de.markers, file = "oligo_pglut_de_markers.csv")

# glutamatergic neurons vs astrocytes
pglut_astro.de.markers <- FindMarkers(data, ident.1 = "23 P Glut", ident.2 = "30 Astro-Epen")
head(pglut_astro.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# LOC103787440     0   5.352007 0.939 0.045         0
# LOC100413948     0   5.000467 0.993 0.100         0
# KSR2             0   4.767807 0.971 0.095         0
# SCN3A            0   4.328221 0.957 0.088         0
# ANO4             0   4.711078 0.939 0.075         0
# PLPPR1           0   4.557623 0.921 0.059         0
range(pglut_astro.de.markers$avg_log2FC)
# [1] -8.625408  6.539726
write.csv(pglut_astro.de.markers, file = "pglut_astro_de_markers.csv")

# astrocytes vs oligodendrocytes
astro_oligo.de.markers <- FindMarkers(data, ident.1 = "30 Astro-Epen", ident.2 = "31 OPC-Oligo")
head(astro_oligo.de.markers)
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# BMPR1B       0   6.470252 0.924 0.036         0
# RFX4         0   6.113532 0.896 0.033         0
# GLI3         0   7.028026 0.866 0.013         0
# SLC4A4       0   6.065063 0.962 0.113         0
# ARHGAP26     0   4.884257 0.935 0.096         0
# GASK1A       0   6.239898 0.865 0.030         0
range(astro_oligo.de.markers$avg_log2FC)
# [1] -5.352012  9.268554
write.csv(astro_oligo.de.markers, file = "astro_oligo_de_markers.csv")


# differential accessibility with tracks----------------------------------------

# plot info for well represented cells
DefaultAssay(data) <- "ATAC"
idents.plot <- c("31 OPC-Oligo", "23 P Glut", "30 Astro-Epen")
data <- SortIdents(data)

# find peaks associated with gene activities
library(presto)

oligo_pglut.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "31 OPC-Oligo",
  ident.2 = "23 P Glut",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(oligo_pglut.da.peaks, file = "oligo_pglut_da_peaks.csv")

astro_oligo.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "30 Astro-Epen",
  ident.2 = "31 OPC-Oligo",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(astro_oligo.da.peaks, file = "astro_oligo_da_peaks.csv")

pglut_astro.da.peaks <- FindMarkers(
  object = data,
  ident.1 = "23 P Glut",
  ident.2 = "30 Astro-Epen",
  test.use = 'wilcox',
  min.pct = 0.1
)
write.csv(pglut_astro.da.peaks, file = "pglut_astro_da_peaks.csv")

# expressed genes linked to peaks
unique(Links(data[['ATAC']])$gene)
# [1] "PRUNE2"  "GLIS3"   "SH3GL2"  "SLC24A2" "ZNRF3"   "FSTL4"   "MAN2A1"  "EDIL3"  
# [9] "VCAN"    "ADCY2"   "GPM6A"   "TMEM144" "BMPR1B"  "SCD5"    "SLC4A4"  "LDB2"   
# [17] "SLC35F1" "MAP7"    "QKI"     "EPN2"    "PITPNC1" "MSI2"    "SV2B"    "CHN1"   
# [25] "R3HDM1"  "PLCL1"   "SATB2"   "AGAP1"   "PIP4K2A" "RNF220"  "TTLL7"   "HECW1"  
# [33] "NXPH1"   "PTPRZ1"  "NELL2"   "TMTC2"   "RFX4"    "MEGF11"  "NOVA1"   "PLEKHH1"
# [41] "SMOC1"   "SLC1A2"  "GRIN2A"  "XYLT1"   "PRKCB"   "SHTN1"   "DLGAP2"  "DOCK5"  
# [49] "AQP4"    "DTNA"    "CDH20"   "PHLPP1"  "MBP"     "DGKG"    "ATP13A4" "PXK"    
# [57] "LRIG1"   "ST18"    "PREX2"   "MMP16"   "KHDRBS3" "ENPP2"   "TRPS1"   "TF"     
# [65] "KAT2B"   "ARPP21"  "ATP1A2"  "TNR"     "CNKSR2" 


# coverage plot for gene recently associated with schizophrenia 
CoveragePlot(
  object = data,
  assay = "ATAC",
  expression.assay = "GENE_ACTIVITIES",
  region = "GRIN2A",
  extend.upstream = 1000,
  extend.downstream = 1000,
  idents = idents.plot,
  annotation = TRUE
)

SaveSeuratRds(object = data, file = "03_multiome_data_filtered.Rds")
