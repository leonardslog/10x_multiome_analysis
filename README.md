# Processing and Analysis of 10X Multiomic data

## Overview

To familiarize myself with RNA and ATAC-seq analysis workflows, I opted to work raw data publicly available via the Allen Institute's <a href="https://brain-map.org/bkp">Brain Knowledge Platform</a>. These analyses were to be carried out on a laptop with 8 cores and 64GB of memory, so to address computational resource limitations, I analyzed 10X Multiomic data generated from brainstem tissue in the common marmoset *Callithrix jacchus*.

Much of the following workflow and code are adapted from the following tutorials/vignettes: \
https://stuartlab.org/signac/articles/overview \
https://satijalab.org/seurat/articles/get_started_v5_new \
https://ngs101.com/tutorials/ \
https://github.com/mousepixels/sanbomics_scripts

## Pre-processing


## Cell Annotation

I annotated cells using the MapMyCells platform with the 10x Whole Mouse Brain Taxonomy (CCN20230722) as the reference.


## Differential Expression analysis


## Differential Accessibility analysis

## References

### Data

Dataset: Fenna M. Krienen, Kirsten M. Levandowski, Heather Zaniewski, Ricardo C.H. del Rosario, Margaret E. Schroeder, Melissa Goldman, Alyssa Lutservitz, Qiangge Zhang, Katelyn X. Li, Victoria F. Beja-Glasser, Jitendra Sharma, Tay Won Shin, Abigail Mauermann, Alec Wysoker, James Nemesh, Seva Kashin, Josselyn Vergara, Gabriele Chelini, Jordane Dimidschstein, Sabina Berretta, Ed Boyden, Steven A. McCarroll, Guoping Feng (2023) A Marmoset Brain Cell Census Reveals Persistent Influence of Developmental Origin on Neurons. [Dataset]. Available from https://assets.nemoarchive.org/dat-1je0mn3

UCSC Genome Browser assembly ID: "calJac4". Perez et al. The UCSC Genome Browser database: 2025 update. Nucleic Acids Research 2025 PMID: 39460617, DOI: 10.1093/nar/gkae974 

### Literature

Ezran, Camille, et al. "Mouse lemur cell atlas informs primate genes, physiology and disease." Nature 644.8075 (2025): 185-196. 

Gontarz, Paul, et al. "Comparison of differential accessibility analysis strategies for ATAC-seq data." Scientific reports 10.1 (2020): 10150.

Krienen, Fenna M., et al. "A marmoset brain cell census reveals regional specialization of cellular identities." Science advances 9.41 (2023): eadk3986. 

Teo, Alan Yue Yang, et al. "Best practices for differential accessibility analysis in single-cell epigenomics." Nature communications 15.1 (2024): 8805. 

Yan, Feng, et al. "From reads to insight: a hitchhiker’s guide to ATAC-seq data analysis." Genome biology 21.1 (2020): 22. 

Yao, Zizhen, et al. "A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain." Nature 624.7991 (2023): 317-332. 

