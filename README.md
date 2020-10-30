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
## Steps Taken to Reach Milestone 1
### Alignment and indexing of FASTQ data with BWA
* The reference file was obtained using the following command: `wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
* The file was unzipped and indexed using BWA: `bwa index -p hg38bwaidx -a bwtsw hg38.fa`
* I then renamed `hr38.fa` as `hg38bwaidx.fa` so all index files have the same name
* The paired reads were aligned to the reference genome using BWA. Sample usage: `bwa mem -M -t 4 hg38bwaidx /scratch/trio/C4RCD_F380_C1_1_1339313_genedx_L001_R1_001.fastq.gz /scratch/trio/C4RCD_F380_C1_1_1339313_genedx_L001_R2_001.fastq.gz > C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sam`
* Three sam files were produced from the previous step, and these were then converted to bam files
* These files were then sorted this step might not be necessary if reproduced since sorting occurs again later
* A read group was added to each file, which is needed for GATK analysis.
* These files were then indexed (it may be possible for this step to be skipped if reproduced, since indexing was done again)
### Mark Duplicates using Picard
* The MarkDuplicates option of Picard was used for this step. Example usage: `java -jar picard/build/libs/picard.jar MarkDuplicates -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.rg.bam -O C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam -M C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.txt`
### Call Variants Using GATK
* Each bam file outputted from the previous step was indexed prior to variant calling.
* Variants were called individually for the three bam files, as well as jointly.
    * Example usage for individual variant calling: `java -jar $GATK HaplotypeCaller -R hg38bwaidx.fa -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam -O C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.indivcall.vcf` 
