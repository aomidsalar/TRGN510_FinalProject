# TRGN 510 Final Project
## Author: Audrey Omidsalar
## Overview of Project
This project will go through a pipeline for analyzing exome sequencing data. The case of interest is a child with neurological disorders characterized by the following features:
* Global development delay - walking by age 5
* Difficulty speaking, though literate
* Hypotonia
* Cleft palate
* Abnormal MRI with reduced white matter
* Tethered spinal cord
* Normal growth, normal spine, above average height

Starting with exome FASTQ files for the child, mother, and father, this pipeline will index and align the files; mark duplicates; call variants, individually and jointly; annotate variants; and finally develop a python script that will filter based on options given, giving a final output of a tab-delimited file containing the filtered results.
## Data
The data for this project includes FASTQ files for the child of interest, father, and mother. This was obtained through the /scratch/trio/ directory in the TRGN server.
## Milestone 1
Completion of the first milestone will require the following steps:
* Alignment and indexing of FASTQ data with BWA
* Marking duplicates using Picard
* Calling variants using GATK
* Annotating variants with SNPEFF
## Milestone 2
The second milestone consists of building a python script which will filter results based on options given from the command-line. The options it will filter for include:
* `--denovo` Showing variants only in the child
* `--impact=LOW,MODERATE,HIGH` Filters for variants which are low, moderate, or high
* `--quality` Providing variants above a given quality score
## Deliverable 
The final product will be a combination of command line scripts through the TRGN server as well as Python scripts from the Jupyter notebook.
