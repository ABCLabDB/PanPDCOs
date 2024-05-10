
# Repository Description

This repository contains the source codes used in the study titled **"The Origin of Patient-Derived Cancer Organoids from Pathologically Undiagnosed Specimen in Patients with Pancreatobiliary Cancers"**, which has been submitted to a journal.

## Scripts

- `240509_scRNA-seq_mainWork`: R code for analyzing PDAC organoid scRNA-Seq data and generating plots for the manuscript.
- `240510_CancerPanel_mainWork.R`: R code for analyzing cancer panel-seq data and generating plots for the manuscript.

## Data

### 1. scRNASeq
- **Metadata.txt**: Contains patient information for PDAC organoid scRNA-Seq data.
- **Zip files**: Each zip file contains the following processed data from the analysis of PDAC organoid scRNA-Seq:
  - **MolsPerCell.csv**: Contains the number of molecules detected per cell. The file is organized in rows, with each row representing a unique cell identifier followed by the total number of mRNA molecules detected in that cell.
  - **ReadsPerCell.csv**: Includes the number of sequencing reads obtained per cell. Each row lists a cell identifier followed by the count of sequencing reads for that cell.

### 2. CancerPanelSeq
- **PheatmapMetadata.txt**: Contains the metadata for cancer panel-seq data.
- **AnnotatedVariants.xlsx**: Contains the variant information for the cancer panel-seq data.

## Citation
TBA
