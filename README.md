# Virus Evolution Lab Avian Influenza Data Pipeline and Analyses

This repository contains notebooks and scripts related to the Nelson (Virus Evolution) Lab's avian influenza research at the National Institutes of Health's National Library of Medicine.

## Purpose

Investigate outbreaks and conduct phylodynamic analyses of avian influenza.

## Authorship
All code in this repository was written by Alexander Maksiaev, with some exceptions:
* [Tree-relabeling script](./metadata-fixing/tree_relabel_keep_bootstrap_aim.py): Alvin Crespo-Bellido
* [Human case-by-genotype lollipop plot](./avian-flu-cats/human_case_by_genotype_lollipop_plot.R): Alvin Crespo-Bellido
* Some [GenoFLU](https://github.com/USDA-VS/GenoFLU) scripts in the archive directory: USDA

## Repository Structure

Avian_Flu/ <br>
│ <br>
├── avian-flu-pipeline/ <br>
│   ├── bash-scripts/ <br>
|   ├── pipeline_documentation.docx <br>
│   ├── Step1_avian_flu_NCBI_Virus_v6.ipynb <br>
│   ├── Step2_avian_flu_Andersen_v6.ipynb <br>
│   ├── Step3_avian_flu_GISAID_v6.ipynb <br>
│   ├── concatenator_v2.ipynb <br>
│   ├── rejects.ipynb <br> 
│   └── utils.py <br>
│ <br>
├── metadata_fixing/ <br>
├── avian-flu-cats/ <br>
├── genotype-a3/ <br>
├── project-specific/ <br>
├── references/ <br>
├── archive/ <br>
└── .gitignore <br>

### avian-flu-pipeline

This pipeline integrates genomic flu sequence data and associated metadata from multiple sources (GenBank, SRA, and GISAID) and prepares structured datasets for downstream phylodynamic analysis. These datasets are updated weekly as an ongoing effort to track the spread of avian influenza since November 2021 within the western hemisphere (Antarctica, North America, and South America), though can be used to create datasets in other locations and time periods as well. 

Specific instructions on how to use this pipeline are detailed in [pipeline_documentation.docx](./avian-flu-pipeline/pipeline_documentation.docx)

### metadata_fixing

One-off scripts designed to update incorrect or out-dated metadata from the avian flu pipeline.

### avian-flu-cats

Scripts related to investigating the spread of avian influenza in felines in the United States, as well as a poster detailing results.

Results: potential transmission sources of avian influenza (specifically genotype B3.13) to cats in the United States were identified: contaminated raw milk from B3.13-infected **dairy cows**, **raw pet food** from B3.13-infected poultry, and caught B3.13-infected **wild birds**.

### genotype-a3

Scripts specific to datasets created for the A3 genotype.

### project-specific

Scripts specific to other projects.

### references

Lookup tables used as references in the avian influenza pipeline, including a reference for which animal species fall into which host type, and abbreviations of administrative regions for various countries (e.g. states within the United States). 

### archive

Deprecated or out-dated workflows. 

