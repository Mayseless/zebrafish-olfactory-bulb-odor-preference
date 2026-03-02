\# mayseless\_et\_al\_2026



This repository contains the analysis scripts and workflows for the study: 

\*\*"A transcriptionally defined interneuron module selectively mediates positive odor valence"\*\*



\## Project Overview

This study establishes a multimodal atlas of the larval zebrafish olfactory bulb (OB) by integrating:

\* \*\*Single-cell RNA sequencing (scRNA-seq):\*\* Resolving transcriptional diversity of the telencephalon.

\* \*\*Whole-mount multiplexed spatial transcriptomics (weMERFISH):\*\* Mapping 68 genes in 3D to resolve spatial subdomains within the OB.

\* \*\*Population-scale activity mapping (CaMPARI):\*\* Linking odor-evoked activity patterns to behavioral valence.

\* \*\*Behavioral Assays:\*\* Quantifying odor preferences using a flow-based two-choice navigation assay.



The integration of these techniques identified a specific dopaminergic, short-axon–like interneuron module (marked by \*th\*) that is essential for positive odor valence (attraction).



\## Repository Structure

The scripts are organized by the experimental modality described in the paper:

\- `1\_Single\_cell/`: Seurat-based clustering and cell-type classification.

\- `2\_Spatial\_transcriptomics/`: Probe design and 3D spatial alignment (Tangram).

\- `3\_Behavior/`: Analysis of flow-based two-choice navigation.

\- `4\_CaMPARI/`: Processing of odor-evoked activity maps.

\- `5\_CaMPARI\_Seq\_\&\_HCR\_Quantification/`: Integration of activity and molecular identity.



\## Methods \& Key Software

As detailed in the supplemental methods, the following tools were used:

\- \*\*scRNA-seq:\*\* CellRanger (vGRCz11) and Seurat (v4.3.0).

\- \*\*Spatial Mapping:\*\* Tangram for projecting scRNA-seq data into 3D space.

\- \*\*Clustering:\*\* Louvain algorithm with resolution optimization via silhouette scores.



\## Citation

If you use these scripts or the data provided, please cite:

> Mayseless, O., Navajas Acedo, J., Wan, Y., Hahaut, V., Picelli, S., Friedrich, R.W., \& Schier, A.F. (2026). A transcriptionally defined interneuron module selectively mediates positive odor valence.



---

\*Note: Large data files (FASTQ, BAM, and CellRanger 'outs') are hosted on \[Zenodo/OSF/Other Link] due to GitHub size constraints.\*

