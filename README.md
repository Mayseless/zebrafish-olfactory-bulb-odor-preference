# Mayseless et al 2026

This repository contains the analysis scripts and workflows for the study:

Mayseless et al. (2026)
**"Odor preference maps to cohesive transcriptional domains in the olfactory bulb"**





## Project Overview

This project establishes a multimodal atlas of the larval zebrafish olfactory bulb (OB), integrating behavioral, functional, and molecular datasets to understand how sensory inputs are transformed into behaviorally relevant outputs.

The study combines:

* Single-cell RNA sequencing (scRNA-seq)

Resolves transcriptional diversity and defines neuronal subtypes.



* Whole-mount spatial transcriptomics

Maps gene expression in 3D to localize transcriptional subtypes within the OB.



* Population-scale activity mapping (CaMPARI)

Captures odor-evoked activity patterns across the entire brain.



* Behavioral assays (two-choice flow paradigm)

Quantifies odor-driven preference and navigation strategies.



* CaMPARI-seq and HCR validation

Links activity-defined neurons to transcriptional identity and spatial localization.





## Repository Structure

The scripts are organized by the experimental modality described in the paper:

* `1_Single_cell/`: Seurat-based clustering and cell-type classification.
* `2_Spatial_transcriptomics/`: Spot detection and cellular assignment, integration with singe cell data (Tangram) and spatial analysis.
* `3_Behavior/`: Analysis of flow-based two-choice navigation.
* `4_CaMPARI/`: Processing of odor-evoked activity maps.
* `5_CaMPARI_Seq_&_HCR_Quantification/`: Integration of activity and molecular identity.





## Notes

This repository accompanies a research manuscript and is intended for transparency and reproducibility, not as a fully packaged software tool.

Code may require adaptation depending on local environment and data structure.



## Citation

If you use these scripts or the data provided, please cite:

> Mayseless, O., Navajas Acedo, J., Wan, Y., Hahaut, V., Picelli, S., Friedrich, R.W., \& Schier, A.F. (2026).Spatial transcriptional modules in the zebrafish olfactory bulb reflect odor preference.

