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
* The reference file was obtained using the following command:

 `wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`
* The file was unzipped and indexed using BWA: 

`bwa index -p hg38bwaidx -a bwtsw hg38.fa`
* I then renamed `hr38.fa` as `hg38bwaidx.fa` so all index files have the same name
* The paired reads were aligned to the reference genome using BWA.
    * `bwa mem -M -t 4 hg38bwaidx /scratch/trio/C4RCD_F380_C1_1_1339313_genedx_L001_R1_001.fastq.gz /scratch/trio/C4RCD_F380_C1_1_1339313_genedx_L001_R2_001.fastq.gz > C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sam`
    * `bwa/bwa mem -M -t 4 hg38bwaidx /scratch/trio/C4RCD_F380_D1_1_1339315_genedx_L001_R1_001.fastq.gz /scratch/trio/C4RCD_F380_D1_1_1339315_genedx_L001_R2_001.fastq.gz > C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sam`
    * `bwa/bwa mem -M -t 4 hg38bwaidx /scratch/trio/C4RCD_F380_M1_2_1339314_genedx_L001_R1_001.fastq.gz /scratch/trio/C4RCD_F380_M1_2_1339314_genedx_L001_R2_001.fastq.gz > C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sam`
* Three sam files were produced from the previous step, and these were then converted to bam files
    * `samtools view -S -b C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sam > C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.bam`
    * `samtools view -S -b C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sam > C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.bam`
    * `samtools view -S -b C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sam > C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.bam`
* These files were then sorted
    * `samtools sort C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.bam -o C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.bam`
    * `samtools sort C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.bam -o C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.bam`
    * `samtools sort C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.bam -o C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.bam`
* A read group was added to each file, which is needed for GATK analysis.
    * `java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.bam -O C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.rg.bam -RGID C1_1_1339313 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM C1_indiv`
    * `java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups -I C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.bam -O C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.rg.bam -RGID D1_1_1339315 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM D1_indiv`
    * `java -jar picard/build/libs/picard.jar AddOrReplaceReadGroups -I C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.bam -O C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.rg.bam -RGID M1_2_1339314 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM M1_indiv`
* These files were then indexed \(it may be possible for this step to be skipped if reproduced, since indexing was done again\)
    * `samtools index C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.rg.bam`
    * `samtools index C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.rg.bam`
    * `samtools index C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.rg.bam`
### Mark Duplicates using Picard
* The MarkDuplicates option of Picard was used for this step.
    * `java -jar picard/build/libs/picard.jar MarkDuplicates -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.rg.bam -O C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam -M C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.txt`
    * `java -jar picard/build/libs/picard.jar MarkDuplicates -I C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.rg.bam -O C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.bam -M C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.txt`
    * `java -jar picard/build/libs/picard.jar MarkDuplicates -I C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.rg.bam -O C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.deduped.bam -M C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.deduped.txt`
### Call Variants Using GATK
* Each bam file outputted from the previous step was indexed prior to variant calling.
    * `samtools index C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam`
    * `samtools index C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.bam`
    * `samtools index C4RCD_F380_M1_2_133C9314_genedx_L001_sequence12_pe.sorted.deduped.bam`
* Variants were called individually for the three bam files, as well as jointly.
    * `java -jar $GATK HaplotypeCaller -R hg38bwaidx.fa -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam -O C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.indivcall.vcf`
    * `java -jar $GATK HaplotypeCaller -R hg38bwaidx.fa -I C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.bam -O C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.indivcall.vcf`
    * `java -jar $GATK HaplotypeCaller -R hg38bwaidx.fa -I C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.deduped.bam -O C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.deduped.indivcall.vcf`
    * Joint Call: `java -jar $GATK HaplotypeCaller -R hg38bwaidx.fa -I C4RCD_F380D_C1_1_1339313_genedx_L001_sequence12_pe.sorted.deduped.bam -I C4RCD_F380_D1_1_1339315_genedx_L001_sequence12_pe.sorted.deduped.bam -I C4RCD_F380_M1_2_1339314_genedx_L001_sequence12_pe.sorted.deduped.bam -O 3samples_jointcall.vcf` 
### Annotating Variants with SnpEff
* SnpEff was used to annotate variants from the `3samples_jointcall.vcf` file generated in the previous step. The reference genome used for SnpEff was the most recent version available to download on the software, GRCh38.99. SnpEff was run in canonical mode.
    * `java -Xmx4G -jar /scratch/aomidsal/TRGN510/FinalProject/snpEff/snpEff.jar -classic -canon -c /scratch/aomidsal/TRGN510/FinalProject/snpEff/snpEff.config GRCh38.99 3samples_jointcall.vcf > 3samples_jointcall.snpEff.vcf`
### Milestone 1 Deliverables
* The deliverable for Milestone 1 consist of the annotated vcf file `canon.3samples_jointcall.snpEff.vcf`. It can be found on the trgn server at the path listed below. There are four files in this directory: the joint call vcf file which was the input for snpEFF \(`canon.3samples_jointcall.vcf`\), the annotated vcf file \(`canon.3samples_jointcall.snpEff.vcf`\), a snpEff genes txt file \(`snpEff_genes.txt`\), and a snpEff summary html file \(`snpEff_summary.html`\)

`/scratch/aomidsal/TRGN510/FinalProject/snpeffcanon`

## Steps Taken to Reach Milestone 2
### Building Python Script
This script was created on JupyterLab, and filters based on the options previously described. It also has an option for an input file, and it prints the tab-delimited output file to stdout so that the user can redirect the output using a file name of his or her choice. 

All of the filtering was done using pandas and regular expressions. Filtering variants based on impact was done by searching for "HIGH", "MODERATE", or "LOW" in the "Info" column of the vcf file. Filtering based on quality is done by searching for variants with a quality score greater than or equal to the integer inputted by the user in the quality option. Filtering for denovo variants searches for those that are only present in the child. In the file used for testing, there were 12 columns in total, with the last three columns containing information about the genotype of the child, dad, and mom, in that order. Therefore, filtering for denovo variants searches for variants that the child possesses that both of the parents do not.

The input file used was the vcf file generated in the previous step `canon.3samples_jointcall.snpEff.vcf.txt`. A `.txt ` was added to the end of the file during the process of transferring the file from the TRGN server onto JupyterLab in order to avoid issues with the computer recognizing the file as a contact list.
### Usage
This python script, `mutationcalling.py`, has 4 options:
* An input file has to be given, using the option `-f` or `--f` followed by the file name
* The option `--denovo` will filter for denovo variants, those where the individual in column 10 is heterozygous or homozygous for the alternate allele, while the individuals in columns 11 and 12 are homozygous for the reference allele. 
* The option `--impact=` can be used to filter by impact of LOW, MODERATE, or HIGH.
* The option `--quality` followed by an integer can be used to filter for variants greater than or equal to the quality score value you entered.
* When running the program, the output tsv file will be written to stdout, so when running the script, redirect the output using `>` and enter your desired file name. 

The output file will be written with columns matching the last line of commented text in the original input file, which will be the column names of your input vcf file. Sample usage:

`!python mutationcalling.py --impact=HIGH,MODERATE --denovo -f canon.3samples_jointcall.snpEff.vcf.txt > filteredvcf.tsv` This will filter for denovo variants that are have a HIGH or MODERATE impact, and the output file will be called `filteredvcf.tsv`
### Determining Mutation
This final project was done with the end goal of finding the denovo mutation that the child possesses which is causing the disease described in the project overview. After running the initial fastq files through the pipeline described in Milestone 1 and filtering using the python script from Milestone 2, I generated two files with my predicted genes of interest. One was generated by filtering for denovo variants that are high impact, and the other was generated by filtering for denovo variants with a quality score greater than or equal to 500. Each of these files had two entries, and after looking up those genes, it was determined that the child has a non-synonymous mutation in the protein coding gene ENST00000400181, otherwise known as KDM1A, and has Kabuki Syndrome. At position 23068566 of chromosome 1, the reference allele is G and the alternate allele is A -- the child is heterozygous, while both parents are homozygous for the reference allele; this is a moderate mutation and was found by filtering the input file for denovo variants with a quality score greater than or equal to 500.

### Dependencies
* BWA
* Samtools
* Picard
* GATK
* snpEff
* Python
    * sys
    * OptionParser
    * Regular Expressions
    * Numpy
    * Pandas
    * CSV
    * FileInput
