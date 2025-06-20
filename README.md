## Analysis pipelines for the pilot phase (V0) of Pig *Cis*-elements Atlas (PCisAtlas)

### 1. Introduction

‘As part of the international Farm animal Genotype-Tissue Expression ([FarmGTEx](https://www.farmgtex.org/)) project, we present Pig *Cis*-elements Atlas ([PCisAtlas](https://PCisAtlas.farmgtex.org/)), a comprehensive multi-tissue developmental atlas of single-cell chromatin accessibility in pigs. This resource encompasses 2.1 million nuclei, profiled from 53 distinct tissues across 37 pigs spanning 12 developmental stages, representing key milestones from early embryogenesis through fetal growth and postnatal development to adulthood. Using PCisAtlas, we systematically annotated 160 distinct cell types and identified 1.77 million candidate cis-regulatory elements (cCREs), significantly expanding the current functional annotation of the pig genome. This unprecedented resource provides a foundational framework for dissecting cell-type-specific regulatory programs and cellular differentiation trajectories. Additionally, PCisAtlas enabled the development of deep learning models to predict the effects of genetic variants across developmental stages and cell types. By integrating these predictions with genome-wide associations (GWAS) of six complex traits in pigs, we prioritized trait-relevant cell types and developmental stages. Furthermore, by examining 444,440 fetal and adult cells in humans, we explored the functional conservation of regulatory elements between humans and pigs. Finally, we leveraged PCisAtlas to interpret GWAS results of 247 complex traits and diseases in humans, providing novel insights into regulatory mechanisms underlying complex traits and diseases at both the cell- and development-specific resolution.’

![PCisAtlas](https://github.com/FarmGTEx/PCisAtlas_pipeline_V0/blob/main/Image/PCisAtlas.png)

---
### 2. Table of Contents in Code/
#### - 01.fastq_process
> Processing of raw FASTQ files, including demultiplexing, read mapping and filtering, and the generation and filtering of fragment/metadata files.
#### - 02.initialize_ArchR
> Basic workflow for initializing the ArchR project, performing cell clustering and visualization, and calling peaks.
#### - 03.NMF analysis
> Identification of cCRE modules using non-negative matrix factorization (NMF).
#### - 04.cell_state_transition_tree
> Construction of the cell state transition tree.
#### - 05.GWAS_enrichment
> Enrichment analysis of pig/human GWAS signals.
#### - 06.pig_human_integration
> Cross-species integration of single-cell data from pig and human.

---
### 3. Citation
**Pending...**
