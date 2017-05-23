# DTU-27626-cancer
Cancer-related exercises for [DTU course 27626](http://www.cbs.dtu.dk/courses/27626/programme.php) 

*WORK IN PROGRESS*
 
These exercises will guide you through all steps starting from raw data (FASTQ files) and resulting in a list of somatic mutations point mutations and copy number changes.  Also, we will perform some analysis (in R) of the resulting data.

Estimated time:  XX hours

## Prerequisites

These exercises are tested with:
* Picard v. 2.9.1
* GATK  v. 3.7
* BWA  v. 0.7.15
* Samtools v. 1.4.1
* Sequenza v. 2.1.2
* R v. 3.4.0


## The data

You will analyze whole-exome sequencing data from a pancreatic tumor and matched normal tissue.

The data used in this exercise has been released for scientific and educational use by the [Texas Cancer Research Biobank](http://txcrb.org/data.html) and is fully described in [this paper](https://www.nature.com/articles/sdata201610).

Please note the Conditions of Data Use:

> By downloading or utilizing any part of this dataset, end users must agree to the following conditions of use:
> * No attempt to identify any specific individual represented by these data or any derivatives of these data will be made.
> * No attempt will be made to compare and/or link this public data set or derivatives in part or in whole to private health information.
> * These data in part or in whole may be freely downloaded, used in analyses and repackaged in databases.
> * Redistribution of any part of these data or any material derived from the data will include a copy of this notice.
> * The data are intended for use as learning and/or research tools only.
> * This data set is not intended for direct profit of anyone who receives it and may not be resold.
> * Users are free to use the data in scientific publications if the providers of the data (Texas Cancer Research Biobank and Baylor College of Medicine Human Genome Sequencing Center) are properly acknowledged.

The raw data files are located on the server at /xxx/xxx/xxx


## Cleanup and alignment (FASTQ file -> BAM file)

1. Step 1

        bwa mem -M -t 12 -R $ReadGoupID $REFF <(bzcat $f1) <(bzcat $f2) | samtools view -b -@ 1 - > $outname

2. Step 2

3. Step 3

## Somatic mutation calling (BAM file -> VCF file)

## Inference of copy number profile

## Data visualization

