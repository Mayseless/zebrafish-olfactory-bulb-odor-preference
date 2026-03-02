# mayseless_et_al_2026

This repository contains the analysis scripts and workflows for the study: 
**"A transcriptionally defined interneuron module selectively mediates positive odor valence"**

## Project Overview
This study establishes a multimodal atlas of the larval zebrafish olfactory bulb (OB) by integrating:
* **Single-cell RNA sequencing (scRNA-seq):** Resolving transcriptional diversity of the telencephalon.
* **Whole-mount spatial transcriptomics:** Mapping genes in 3D space to resolve spatial subdomains within the OB.
* **Population-scale activity mapping (CaMPARI):** Linking odor-evoked activity patterns to behavioral valence.
* **Behavioral Assays:** Quantifying odor preferences using a flow-based two-choice navigation assay.

The integration of these techniques identified a specific dopaminergic, short-axon–like interneuron module (marked by *th*) that is essential for positive odor valence (attraction).

## Repository Structure
The scripts are organized by the experimental modality described in the paper:
- `1_Single_cell/`: Seurat-based clustering and cell-type classification.
- `2_Spatial_transcriptomics/`: Spot detection and cellualr assignment, integration with singe cell data (Tangram) and spatial analysis.
- `3_Behavior/`: Analysis of flow-based two-choice navigation.
- `4_CaMPARI/`: Processing of odor-evoked activity maps.
- `5_CaMPARI_Seq_&_HCR_Quantification/`: Integration of activity and molecular identity.


## Citation
If you use these scripts or the data provided, please cite:
> Mayseless, O., Navajas Acedo, J., Wan, Y., Hahaut, V., Picelli, S., Friedrich, R.W., & Schier, A.F. (2026). A transcriptionally defined interneuron module selectively mediates positive odor valence.
