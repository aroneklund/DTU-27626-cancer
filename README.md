# DTU-27626-cancer
Cancer-related exercises for [DTU course 27626](http://www.cbs.dtu.dk/courses/27626/programme.php) 

Marcin Krzystanek and Aron Eklund
 
These exercises will guide you through all steps starting from raw data (FASTQ files) and resulting in a list of somatic mutations point mutations and copy number changes.  Also, we will perform some analysis (in R) of the resulting data.

Estimated time:  2 hours

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

The raw data files are located on the server at /home/27626/exercises/cancer

Important: since data is vast and our resources are limited we will not run the alignment and full mutation calling. We do provide the code needed 

### 1.1 Take a first look at the data

        ls /home/27626/exercises/cancer
        bzcat /home/27626/exercises/cancer/TCRBOA2-N-WEX.read2.fastq.bz2 | head

Q1: Is your data single or paired end? What type would you prefer for cancer DNA sequencing and why?

## 2. Data Pre-processing

        ### Define bash variables for brevity:
        f1n=/home/27626/exercises/cancer/TCRBOA2-N-WEX.read1.fastq.bz2
        f2n=/home/27626/exercises/cancer/TCRBOA2-N-WEX.read2.fastq.bz2
        f1t=/home/27626/exercises/cancer/TCRBOA2-T-WEX.read1.fastq.bz2
        f2t=/home/27626/exercises/cancer/TCRBOA2-T-WEX.read2.fastq.bz2
        HREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set
        FREFF=/home/27626/exercises/cancer/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fa 
        IREFF=/home/27626/exercises/cancer/human_GRCh38/Indel_refs/mills_gold.b38.vcf
        SREFF=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/1000G.snps.b38.vcf
        cosmicREFF=/home/27626/exercises/cancer/human_GRCh38/cosmic/CosmicCodingMuts_chr_sorted.vcf
        GATK=/home/27626/exercises/cancer/programs/GenomeAnalysisTK.jar
        PICARD=/home/27626/exercises/cancer/progrpicard-2.jar
        TRIM_GALORE=/home/27626/exercises/cancer/programs/trim_galore
        outdir=`pwd`


### 2.1 [DO NOT RUN] Read quality trimming and FastQC report using [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

Next step in which you are going computations is point 3 but we left the code needed for preparation of the data below so that you can use it later in the course project or in your own work should you work with cancer patient DNA sequencing data.

        ### Arguments to be passed to FastQC
        args="'--outdir ${outdir}'"
        
        ### Trim reads with trim_galore wrapper, produce both fastqc and trimming reports
        $TRIM_GALORE --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1n $f2n
        $TRIM_GALORE --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1t $f2t

Q2: What does --quality 20 argument mean? Get help by running:
        
        $TRIM_GALORE --help

Set up new variables for the newly created files. I assume the validated and filtered files were created in your working directory (for me this is /home/27626/exercises/cancer/ so you can find these files there if you need them).

        f1n_val=TCRBOA2-N-WEX.read1.fastq.bz2_val_1.fq.gz
        f2n_val=TCRBOA2-N-WEX.read2.fastq.bz2_val_2.fq.gz
        f1t_val=TCRBOA2-T-WEX.read1.fastq.bz2_val_1.fq.gz
        f2t_val=TCRBOA2-T-WEX.read2.fastq.bz2_val_2.fq.gz


### 2.2 [DO NOT RUN] Alignment and preprocessing before mutation calling.


2.2.1 Step 1 - (PLEASE DO NOT RUN)

Use bwa mem for aligning reads. Tumor sample and normal separately.

Importantly, a Read Group ID line (@RG line) needst to be defined by the user. Mutect2 and other programs in the pipeline below depend on information in this line. Here is one way of constructing it. Optionally, it can have more information then provided below. Please see the [SAM format specification](http://www.samformat.info) if you want to know more.

        ### @RG ID # read group ID, needs to be unique for fastq file due to downstream processing, takes\
        preferrence when used by some programs
        ### @RG SM # sample ID, unique for each tumor and normal sample, not to be confused with patient ID
        ### @RG PL # platform name
        ### @RG LB # library name
        ### @RG PU # Platform unit, needs to be unique for fastq file due to downstream processing, takes\
        preferrence when used by some programs
        ### Let's create an @RG line that we will use when runnig bwa mem alinment
        ReadGoupID_N="\"@RG\tID:TCRBOA2-N-WEX\tSM:TCRBOA2-N-WEX\tPL:ILLUMINA\tLB:libN\tPU:TCRBOA2-N-WEX"\"
        ReadGoupID_T="\"@RG\tID:TCRBOA2-T-WEX\tSM:TCRBOA2-T-WEX\tPL:ILLUMINA\tLB:libT\tPU:TCRBOA2-T-WEX"\"

        ### Run bwa mem
        bwa mem -M -t 4 -R $ReadGoupID_N $HREFF $f1n_val $f2n_val \
            | samtools view -Sb -@ 1 - > patient2_n.bam 
        bwa mem -M -t 4 -R $ReadGoupID_T $HREFF $f2t_val $f2t_val \
            | samtools view -Sb -@ 1 - > patient2_t.bam

2.2.2 Step 2 - Sort bam files. 

        samtools sort -@ 3 patient2_n.bam -o patient2_n.sorted.bam
        samtools sort -@ 3 patient2_t.bam -o patient2_t.sorted.bam
 
2.2.3 Step 3 - Mark duplicates with picard tools MarkDuplicates. 
Mark PCR duplicates so that they will not introduce false positives and bias in the subsequent analysis.

        mkdir tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=patient2_n.sorted.bam OUTPUT=patient2_n.sorted.dedup.bam METRICS_FILE=patient2_n.metrics.txt \
            TMP_DIR=./tmp
        java -Xmx10G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar $PICARD MarkDuplicates\
            INPUT=patient2_t.sorted.bam OUTPUT=patient2_t.sorted.dedup.bam METRICS_FILE=patient2_t.metrics.txt \
            TMP_DIR=./tmp
        

2.2.4 Step 4 - Index bam.

        samtools index patient2_n.sorted.dedup.bam
        samtools index patient2_t.sorted.dedup.bam

2.2.5 Step 5 - BaseRecalibrator - Part 1
Recalibrate base qualities. Each base in each sequence read comes out of sequencer with certain quality score. Depending on machine used for sequrencing these scores are subjected to various sources of systematic technical error. Base quality score recalibration (BQSR) works by applying machine learning to model these errors empirically and adjust the quality scores accordingly. Here is more information on [BSQR](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php) from authors of the software. 

        ### Run Base Recalibrator 1 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient2_n.recal.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient2_t.recal.table
      
        
2.2.6 Step 6 - BaseRecalibrator - Part 2.

        ### Run Base Recalibrator 2 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient2_n.recal.table -o patient2_n.post_recal_data.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T BaseRecalibrator -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient2_t.recal.table -o patient2_t.post_recal_data.table
        
        
2.2.7 Step 7 - PrintReads.
Recalibrated reads are collected in a new bam file. After this step, the resulting bam file is ready to be processed with MuTect2 - mutation calling program.

        ### Run Base Recalibrator - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK\
            -T PrintReads -nct 4 -R $FREFF -I patient2_n.sorted.dedup.bam -BQSR patient2_n.recal.table \
            -o patient2_n.final.bam
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar $GATK \
            -T PrintReads -nct 4 -R $FREFF -I patient2_t.sorted.dedup.bam -BQSR patient2_t.recal.table \
            -o patient2_t.final.bam


## 3. Somatic mutation calling (BAM file -> VCF file)
Since we do not have time and capacity to run a whole sample during our exercises we will call somatic mutations at chromosome 1 from 50.000.000th to 52.000.000th base pair.
MuTect2 is a somatic mutation caller developed by Broad Institute. MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect [(Cibulskis et al., 2013)](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html) with the assembly-based machinery of [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) which you have already used in the [genotyping exercise](http://www.cbs.dtu.dk/courses/27626/programme.php). The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller but MuTect2 allows varying allelic fraction for each variant. This is necessary because tumors often are hetergeneous (multiclonal), have lower cellularity (purity) than 100%, show gains and losses of parts of the genome. To learn more about Mutect2 follow this link [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)
Mutect2 is computationally intensive so it is recommended to parallelize if possible. One way to achieve it is to split processes by chromosomes.

        ### Set chromosome and location:
        CHR_LOC=chr1:50000000-52000000
        ### Use pre-processed bam files:
        fbn=/home/27626/exercises/cancer/patient2_n.final.bam
        fbt=/home/27626/exercises/cancer/patient2_t.final.bam
        ### Run Mutect2
        time java -Xmx4G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=1 -jar $GATK -T MuTect2 \
            -R $FREFF --dbsnp $SREFF --cosmic $cosmicREFF -I:tumor $fbt \
            -I:normal $fbn -o patient2_t.${CHR_LOC}.mutect2.vcf -L $CHR_LOC
        ### To run a whole genome simply do not use the -L option.

Take a look at the VCF file. Unlike HoaplotypeCaller MuTect2 applies a range of filters to each call by default. For a start try to filter mutational calls by selecting those with MuTect "PASS" annotation.

        cat patient2_t.${CHR_LOC}.mutect2.vcf | grep PASS

You should see this line:
"chr1	50973993	rs746646631	C	T	.	PASS	DB;ECNT=1;HCNT=2;MAX_ED=.;MIN_ED=.;NLOD=33.99;TLOD=7.29	GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1	0/1:129,6:0.044:3:3:0.500:3973,169:66:63	0/0:132,0:0.00:0:0:.:4093,0:75:57"

Explanation of each part of the line above is in the header of the file (use less command to look at it). Importantly, the column starting with 0/0 refers to the normal sample while the one beginning with 0/1 refers to the tumor. After genotype (GT) we have allelic depth (AD) hich is "129,6" (i.e. 126 and 6 for the reference and mutant allele respectively). Then comes allelic frequency which is a fraction of the mutant allele out of all aligned bases in this position.

Q3: Try to search [dbSNP](https://www.ncbi.nlm.nih.gov/snp) for rs746646631. What gene does it belong to? Is this mutation protein-chainging?

Go to [cBIO](http://www.cbioportal.org) portal that contains a collection of large cancer datasets and type the name of the gene that was hit by this mutation in the "Enter Gene Set:" box in the bottom of the page. Press submit. How often is this gene mutated in various cancer types?  

## Inference of tissue of origin

        ### .....

## Inference of copy number profile

        ### .....

## Data visualization

        ### .....

