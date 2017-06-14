# DTU-27626-cancer
Cancer-related exercises for [DTU course 27626](http://www.cbs.dtu.dk/courses/27626/programme.php) 

Marcin Krzystanek and Aron Eklund

These exercises will guide you through all steps starting from raw data (FASTQ files)
and resulting in a list of somatic point mutations and a copy number profile.
Also, we will perform some analysis (in R) of the resulting data.

Estimated time:  2 hours

## Prerequisites

These exercises are tested with:
* Picard v. 2.9.1
* GATK  v. 3.7
* BWA  v. 0.7.15
* Samtools v. 1.4.1
* Sequenza v. 2.1.2
* R v. 3.4.0


## About the data

You will analyze whole-exome sequencing data from a pancreatic tumor and matched normal tissue.

The data used in this exercise has been released for scientific and educational use by the
[Texas Cancer Research Biobank](http://txcrb.org/data.html) and is fully described in 
[this paper](https://www.nature.com/articles/sdata201610).

Please note the Conditions of Data Use:

> By downloading or utilizing any part of this dataset, end users must agree to the following conditions of use:
> * No attempt to identify any specific individual represented by these data or any derivatives of these data will be made.
> * No attempt will be made to compare and/or link this public data set or derivatives in part or in whole to private health information.
> * These data in part or in whole may be freely downloaded, used in analyses and repackaged in databases.
> * Redistribution of any part of these data or any material derived from the data will include a copy of this notice.
> * The data are intended for use as learning and/or research tools only.
> * This data set is not intended for direct profit of anyone who receives it and may not be resold.
> * Users are free to use the data in scientific publications if the providers of the data (Texas Cancer Research Biobank and Baylor College of Medicine Human Genome Sequencing Center) are properly acknowledged.

The raw data files are located on the server at `/home/27626/exercises/cancer`

## Somatic point mutation exercise

**IMPORTANT IMPORTANT IMPORTANT** - Since the full procedure takes a long time, 
we will **not** ask you to perform the full alignment and full mutation calling. 
However, for reference, we provide the code needed for the full analysis. 
Thus, you can use this code later in the course project or in your own work, 
should you work with cancer patient DNA sequencing data.

The parts where you should actually run the code include: 1.1, 1.2, 3, and 4


### PART 1. Raw reads: inspection, QC, cleanup

#### 1.1 - Take a first look at the data

        ls /home/27626/exercises/cancer
        bzcat /home/27626/exercises/cancer/TCRBOA2-N-WEX.read2.fastq.bz2 | head

Q1: How long are the reads? Is your data single or paired end? 
What type would you prefer for cancer DNA sequencing, and why?


#### 1.2 - Define some bash variables

        ### Define bash variables for brevity:
        f1n=/home/27626/exercises/cancer/TCRBOA2-N-WEX.read1.fastq.bz2
        f2n=/home/27626/exercises/cancer/TCRBOA2-N-WEX.read2.fastq.bz2
        f1t=/home/27626/exercises/cancer/TCRBOA2-T-WEX.read1.fastq.bz2
        f2t=/home/27626/exercises/cancer/TCRBOA2-T-WEX.read2.fastq.bz2
        HREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set
        FREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fa 
        IREFF=/home/27626/exercises/cancer/human_GRCh38/Indel_refs/mills_gold.b38.vcf
        SREFF=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/1000G.snps.b38.vcf
        dbsnp_ALL=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/dbsnp/All_20160527_chr.vcf
        cosmicREFF=/home/27626/exercises/cancer/human_GRCh38/cosmic/CosmicCodingMuts_chr_sorted.vcf
        GATK=/home/27626/exercises/cancer/programs/GenomeAnalysisTK.jar
        PICARD=/home/27626/exercises/cancer/programs/picard-2.jar
        TRIM_GALORE=/home/27626/exercises/cancer/programs/trim_galore
        outdir=`pwd`


#### 1.3 - Read quality trimming and FastQC report (DO NOT RUN)

We do this using [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)


        ### Arguments to be passed to FastQC
        args="'--outdir ${outdir}'"
        
        ### Trim reads with trim_galore wrapper, produce both fastqc and trimming reports
        $TRIM_GALORE --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1n $f2n
        $TRIM_GALORE --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1t $f2t

Q2: What does the argument `--quality 20` mean? Get help by running:
        
        $TRIM_GALORE --help

Set up new variables for the newly created files. I assume the validated and filtered
files were created in your working directory (for me this is /home/27626/exercises/cancer/
so you can find these files there if you need them).

        f1n_val=TCRBOA2-N-WEX.read1.fastq.bz2_val_1.fq.gz
        f2n_val=TCRBOA2-N-WEX.read2.fastq.bz2_val_2.fq.gz
        f1t_val=TCRBOA2-T-WEX.read1.fastq.bz2_val_1.fq.gz
        f2t_val=TCRBOA2-T-WEX.read2.fastq.bz2_val_2.fq.gz


### PART 2. Alignment and additional preprocessing (DO NOT RUN)

#### 2.1 - Alignment (DO NOT RUN)

We use [bwa mem](https://github.com/lh3/bwa) for aligning reads to the genome. 
We align the tumor sample and normal sample separately.

Importantly, a Read Group ID line (@RG line) must be defined by the user, because Mutect2
and other programs in the pipeline below depend on information in this line. Here we
demonstrate one way of adding the @RG line to the resulting BAM file:

        ### @RG ID # read group ID, needs to be unique for fastq file due to downstream processing, takes\
        preference when used by some programs
        ### @RG SM # sample ID, unique for each tumor and normal sample, not to be confused with patient ID
        ### @RG PL # platform name
        ### @RG LB # library name
        ### @RG PU # Platform unit, needs to be unique for fastq file due to downstream processing, takes\
        preference when used by some programs
        ### Let's create an @RG line that we will use when running bwa mem alignment
        ReadGoupID_N="\"@RG\tID:TCRBOA2-N-WEX\tSM:TCRBOA2-N-WEX\tPL:ILLUMINA\tLB:libN\tPU:TCRBOA2-N-WEX"\"
        ReadGoupID_T="\"@RG\tID:TCRBOA2-T-WEX\tSM:TCRBOA2-T-WEX\tPL:ILLUMINA\tLB:libT\tPU:TCRBOA2-T-WEX"\"

        ### Run bwa mem
        bwa mem -M -t 4 -R $ReadGoupID_N $HREFF $f1n_val $f2n_val \
            | samtools view -Sb -@ 1 - > patient2_n.bam 
        bwa mem -M -t 4 -R $ReadGoupID_T $HREFF $f2t_val $f2t_val \
            | samtools view -Sb -@ 1 - > patient2_t.bam

Optionally, the @RG line can provide additional information; please see the 
[SAM format specification](http://www.samformat.info) as well as [samtools webpage](http://samtools.sourceforge.net) if you want to know more.

#### 2.2 - Sort BAM files (DO NOT RUN)

        samtools sort -@ 3 patient2_n.bam -o patient2_n.sorted.bam
        samtools sort -@ 3 patient2_t.bam -o patient2_t.sorted.bam


#### 2.3 - Mark duplicates (DO NOT RUN)

We use [Picard](https://broadinstitute.github.io/picard/) to mark PCR duplicates so that
they will not introduce false positives and bias in the subsequent analysis.

        mkdir tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=patient2_n.sorted.bam OUTPUT=patient2_n.sorted.dedup.bam METRICS_FILE=patient2_n.metrics.txt \
            TMP_DIR=./tmp
        java -Xmx10G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar $PICARD MarkDuplicates\
            INPUT=patient2_t.sorted.bam OUTPUT=patient2_t.sorted.dedup.bam METRICS_FILE=patient2_t.metrics.txt \
            TMP_DIR=./tmp
        

#### 2.4 - Index the BAM files (DO NOT RUN)

        samtools index patient2_n.sorted.dedup.bam
        samtools index patient2_t.sorted.dedup.bam


#### 2.5 - BaseRecalibrator - Part 1 (DO NOT RUN)

We use [GATK](https://software.broadinstitute.org/gatk/) to recalibrate base quality scores. 

Each base in each sequencing read comes out of the sequencer with an individual quality score. 
Depending on the machine used for sequencing, these scores are subject to various
sources of systematic technical error. Base quality score recalibration (BQSR) works by
applying machine learning to model these errors empirically and adjust the quality scores
accordingly. 

There is more information on BSQR [here](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php). 

        ### Run Base Recalibrator 1 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient2_n.recal.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient2_t.recal.table
      
        
#### 2.6 - BaseRecalibrator - Part 2 (DO NOT RUN)

        ### Run Base Recalibrator 2 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient2_n.recal.table -o patient2_n.post_recal_data.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient2_t.recal.table -o patient2_t.post_recal_data.table
        
        
#### 2.7 - PrintReads (DO NOT RUN)

The recalibrated reads are collected in a new BAM file. 

        ### Run Base Recalibrator - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK\
            -T PrintReads -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -BQSR patient2_n.recal.table \
            -o patient2_n.final.bam
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T PrintReads -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -BQSR patient2_t.recal.table \
            -o patient2_t.final.bam

Now, the resulting BAM files are ready to be processed with MuTect2.


### PART 3. Somatic mutation calling (BAM file -> VCF file)

#### 3.1 - MuTect2

We use [MuTect2][MuTect2], a somatic mutation caller that identifies both SNV and indels.

[MuTect2]: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php

Mutect2 is computationally intensive so we recommend to parallelize if possible. 
One way to achieve this is to split processes by chromosomes.

Since we do not have the time and capacity to process the entire genome during our
exercises, we will call somatic mutations on a small part of chromosome 1, 
from the 50.000.000th to the 52.000.000th base pair.

        ### Set chromosome and location:
        CHR_LOC=chr1:50000000-52000000
        ### Use pre-processed bam files:
        fbn=/home/27626/exercises/cancer/patient2_n.final.bam
        fbt=/home/27626/exercises/cancer/patient2_t.final.bam
        ### Run Mutect2
        java -Xmx4G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -jar $GATK -T MuTect2 \
            -R $FREFF --dbsnp $dbsnp_ALL --cosmic $cosmicREFF -I:tumor $fbt \
            -I:normal $fbn -o patient2_t.${CHR_LOC}.mutect2.vcf -L $CHR_LOC
        ### To process the whole genome, simply omit the -L option.

Take a look at the resulting VCF file. Unlike HaplotypeCaller, MuTect2 applies a
range of filters to each call by default. 

#### 3.2 - Filter the VCF output

For a start, try to filter mutational calls by selecting those with MuTect "PASS" annotation.

        cat patient2_t.${CHR_LOC}.mutect2.vcf | grep PASS

You should see this line:
"chr1	50973993	rs746646631	C	T	.	PASS	DB;ECNT=1;HCNT=2;MAX_ED=.;MIN_ED=.;NLOD=33.99;TLOD=7.29	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1	0/1:129,6:0.044:3:3:0.500:3973,169:66:63	0/0:132,0:0.00:0:0:.:4093,0:75:57"

A brief explanation of each part of the line above is in the header of the VCF file 
(use "less" to look at it). Importantly, the column starting with 0/0 refers to the
normal sample, whereas the column beginning with 0/1 refers to the tumor. 
After genotype (GT) we have allelic depth (AD) which is "129,6" (i.e. 129 and 6 for
the reference and mutant allele respectively). Then comes allelic frequency, which is
a fraction of the mutant allele out of all aligned bases in this position. 
For more information about the MuTect2 output go to [MuTect2][MuTect2].


### PART 4. Interpretation of the resulting somatic mutations

#### 4.1 - Interpretation of somatic mutations

Q3: Try to search [dbSNP](https://www.ncbi.nlm.nih.gov/snp) for rs746646631. 
What gene does it belong to? Is this mutation protein-changing?

Go to [cBioPortal](http://www.cbioportal.org), a website that provides tools to analyze
several large cancer sequencing datasets. Type the name of the gene that was hit by this
mutation in the "Enter Gene Set:" box in the bottom of the page and press submit. 
How often is this gene mutated in various cancer types?  


#### 4.2 - More interpretation of somatic mutations

So far you have processed and analyzed only a small section of chromosome 1.

Now, let's analyze a bigger piece of the genome. Pick your favorite chromosome and
find the corresponding VCF file on the server.  For example, if you choose
chromosome 7, you would use this file:

        ls /home/27626/exercises/cancer/patient2.chr7.mutect2.vcf

Hint: your results will be more interesting if you pick chromosome 
6, 13, 15, 17, 19, 20, 22, or X!

Filter the VCF to retain only the lines marked as "PASS".  

        grep "PASS" /home/27626/exercises/cancer/patient2.chr7.mutect2.vcf > filtered.chr7.vcf

Download the *filtered* VCF to your own computer and submit it to the 
[VEP website](http://www.ensembl.org/Tools/VEP) using default settings. 
When the results become available, look in the "Somatic status" column. Are there
any known cancer mutations?
If you find a known cancer mutation, find its COSMIC identifier 
(COSM######, e.g. COSM4597270) in the "existing variant" column.
Search for your COSMIC identifier in the
[COSMIC database](http://cancer.sanger.ac.uk/cosmic).
In which tissues is this mutation found?


#### 4.3 - Inference of tissue of origin

Next we'll do some analysis on a VCF file containing somatic mutations found throughout
the entire genome:

        ls /home/27626/exercises/cancer/patient2_t.mutect2.vcf

Unlike VEP, TumorTracer requires VCF files to have the header information.
Thus, we will filter this VCF file to retain: 1) header lines (which begin with "#"),
and 2) data lines with a PASS call.
        
        grep -E "^#|PASS" /home/27626/exercises/cancer/patient2_t.mutect2.vcf > filtered.patient2_t.mutect2.vcf

Submit the filtered VCF to the
[TumorTracer server](http://www.cbs.dtu.dk/services/TumorTracer/).
Make sure to specify that this VCF was generated using GRCh38 coordinates.

What tissue does TumorTracer predict?  Is it a confident prediction?


## PART 5. Inference of copy number profile

Sequenza is an R package that can be used to infer the copy number profile of a tumor
specimen. The input is (exome or whole-genome) NGS data from a tumor specimen and a
matched normal specimen.

Links to the Sequenza publication, source code, R package, help forums, etc. can be
found on the [Sequenza home page](http://www.cbs.dtu.dk/biotools/sequenza/).

Sequenza has already been installed on the course server, so you first need to start R, 
and then load the package:

        R
        library(sequenza)

Normally, several steps are required to prepare the data. The raw reads must be aligned
to the genome, and then the BAM file must be processed by a python script (included in
the R package) that extracts the relevant information into a format suitable for analyis
in R. This takes more time and computation power than is reasonable for an exercise,
so we will skip ahead to the fun part.

If you need to preprocess your own data, please see the instructions in the Sequenza
manual, which you can obtain by typing (in R):
        
        vignette("sequenza")
        
####  5.1 Run sequenza

We are going to work on a lung adenocarcinoma specimen (WGS), derived from a
55-year-old female donor, who is still alive and symptom free after a successful
resection. Let's start by loading the data into R:

        load("/home/27626/exercises/cancer/LUADsample_seqz.RData")

This object (seqzDF) is the output of the `sequenza.extract()` function.
Feel free to check its contents (`str(seqzDF)`)! Now we can fit a model to this data.
It will probably keep your computer busy for a few minutes.

        seqz.fit   <- sequenza.fit(seqzDF)

Finally, we create the output files. This creates a subdirectory in your working folder 
called "sequenza-luad" and fills it with 13 files. 

        sequenza.results(seqzDF, seqz.fit, out.dir = "sequenza-luad", sample.id = "luad55")

####  5.2 Interpretation of sequenza results

Check out the files in the "sequenza-luad" directory. (If you are not sure where
this directory is, type `getwd()` in R). You can view the PDFs directly from the server
using `evince`, but it often works better to download the files to your local machine
and view them there.

Questions:

**Q1:** Does the model fit seem reasonable? How do the "alternative" fits look, visually?

**Q2:** Where in the genome do you see the largest copy number?

**Q3:** Do you see any evidence of loss of heterozygosity? Where? (give an example or two)

