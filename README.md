# Zebrafish olfactory bulb odor preference

This repository contains analysis scripts and workflows associated with the manuscript:

**"Odor preference maps to cohesive transcriptional domains in the olfactory bulb"**

## Project overview

This project integrates behavioral, functional, single-cell, and spatial transcriptomic datasets to examine how odorant-evoked activity in the larval zebrafish olfactory bulb relates to odor preference and transcriptionally defined neuronal domains.

The study combines:

- **Single-cell RNA sequencing (scRNA-seq)**  
  Defines transcriptional diversity and neuronal subtypes.

- **Whole-mount spatial transcriptomics**  
  Maps gene expression in 3D and localizes transcriptional subtypes within the olfactory bulb.

- **Population-scale activity mapping (CaMPARI2)**  
  Captures odorant-evoked activity patterns in freely behaving larvae.

- **Behavioral assays**  
  Quantifies odorant-driven preference and navigation strategies using a two-choice flow paradigm.

- **CaMPARI2-seq and HCR validation**  
  Links odorant-activated neurons to transcriptional identity and spatial localization.

## Repository structure

Scripts are organized by experimental modality:

- `1_Single_cell/`  
  Seurat-based clustering, neuronal subtype classification, and marker analysis.

- `2_Spatial_transcriptomics/`  
  Spatial transcriptomic processing, cellular assignment, Tangram-based integration with single-cell data, and spatial organization analyses.

- `3_Behavior/`  
  Analysis of flow-based two-choice behavioral assays and odorant preference.

- `4_CaMPARI/`  
  Processing and analysis of CaMPARI2-based odorant-evoked activity maps.

- `5_CaMPARI_Seq_&_HCR_Quantification/`  
  Analysis of CaMPARI2-seq, HCR validation, and links between activity, molecular identity, and spatial localization.

## Notes

This repository accompanies a research manuscript and is intended to support transparency and reproducibility. It is not a fully packaged software tool.

Scripts may require adaptation depending on the local computing environment, file paths, and data organization.

## Citation

If you use these scripts or associated data, please cite the manuscript once available:

> Mayseless, O., Navajas Acedo, J., Wan, Y., Hahaut, V., Picelli, S., Friedrich, R.W., & Schier, A.F. Odor preference maps to cohesive transcriptional domains in the olfactory bulb.

A permanent citation and DOI will be added upon publication.
