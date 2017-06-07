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
        IREFF=/home/27626/exercises/cancer/human_GRCh38/Indel_refs/mills_gold.b38.vcf
        SREFF=/home/27626/exercises/cancer/human_GRCh38/SNP_refs/1000G.snps.b38.vcf
        cosmicREFF=/home/27626/exercises/cancer/human_GRCh38/cosmic/CosmicCodingMuts_chr_sorted.vcf
        GATK=/home/27626/exercises/cancer/programs/GenomeAnalysisTK.jar
        PICARD=/home/27626/exercises/cancer/progrpicard-2.jar
        TRIM_GALORE=/home/27626/exercises/cancer/programs/trim_galore
        outdir=`pwd`


### 2.1 Read quality trimming and FastQC report using [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

        ### Arguments to be passed to FastQC
        args="'--outdir ${outdir}'"
        
        ### Trim reads with trim_galore wrapper, produce both fastqc and trimming reports
        trim_galore --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1n $f2n
        trim_galore --fastqc --fastqc_args $args --gzip --quality 20 --trim-n --length 50\
        --trim1 --output_dir $outdir --paired $f1n $f2n



### 2.1 Cleanup and alignment (FASTQ file -> BAM file)


1. Step 1 - (PLEASE DO NOT RUN)

Use bwa mem for aligning reads. Tumor sample and normal separately.
 

Importantly and Read Group ID line (@RG line) needst to be defined by the user. Mutect2 and other programs in the pipeline below depend on information in this line. Here is one way of constructing it. Optionally, it can have more information then provided below. Please see the [SAM format specification](http://www.samformat.info) if you want to know more.

        ### @RG ID # read group ID, needs to be unique for fastq file due to downstream processing, takes preferrence when used by some programs
        ### @RG SM # sample ID
        ### @RG PL # platform name
        ### @RG LB # library name
        ### @RG PU # Platform unit, needs to be unique for fastq file due to downstream processing, takes preferrence when used by some programs
        ### Let's create an @RG line that we will use when runnig bwa mem alinment
        ReadGoupID_N="\"@RG\tID:TCRBOA2-N-WEX\tSM:TCRBOA2-N-WEX\tPL:ILLUMINA\tLB:libN\tPU:TCRBOA2-N-WEX"\"
        ReadGoupID_T="\"@RG\tID:TCRBOA2-T-WEX\tSM:TCRBOA2-T-WEX\tPL:ILLUMINA\tLB:libT\tPU:TCRBOA2-T-WEX"\"

        ### Run bwa mem
        bwa mem -M -t 4 -R $ReadGoupID_N $HREFF <(bzcat $f1n) <(bzcat $f2n) \
            | samtools view -Sb -@ 1 - > patient3_n.bam 
        bwa mem -M -t 4 -R $ReadGoupID_T $HREFF <(bzcat $f2t) <(bzcat $f2t) \
            | samtools view -Sb -@ 1 - > patient3_t.bam

2. Step 2 - Sort bam files. 

        samtools sort -@ 3 patient3_n.bam -o patient3_n.sorted.bam
        samtools sort -@ 3 patient3_t.bam -o patient3_t.sorted.bam
 
3. Step 3 - Mark duplicates with picard tools MarkDuplicates. 
Mark PCR duplicates so that they will not introduce false positives and bias in the subsequent analysis.

        
        mkdir tmp
        java -Xmx5G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=6 -jar $PICARD MarkDuplicates\
            INPUT=patient3_n.sorted.bam OUTPUT=patient3_n.sorted.dedup.bam METRICS_FILE=patient3_n.metrics.txt \
            TMP_DIR=./tmp
        java -Xmx10G -Xms1024M  -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar $PICARD MarkDuplicates\
            INPUT=patient3_t.sorted.bam OUTPUT=patient3_t.sorted.dedup.bam METRICS_FILE=patient3_t.metrics.txt \
            TMP_DIR=./tmp
        

4. Step 4 - Index bam.

        samtools index patient3_n.sorted.dedup.bam
        samtools index patient3_t.sorted.dedup.bam

5. Step 5 - BaseRecalibrator - Part 1
Recalibrate base qualities. Each base comes out of sequencer with certain quality call. These qualities are readjusted by BaseRecalibrator after [TO BE FINISHED]

        ### Run Base Recalibrator 1 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATK}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient3_n.recal.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATK}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient3_t.recal.table
      
        
6. Step 6 - BaseRecalibrator - Part 2.

        ### Run Base Recalibrator 2 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATK}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient3_n.recal.table -o patient3_n.post_recal_data.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -BQSR patient3_t.recal.table -o patient3_t.post_recal_data.table
        
        
6. Step 6 - PrintReads.
Recalibrated reads are collected in a new bam file. The new bam file is now ready to be processed with MuTect2 - mutation calling program.

        ### Run Base Recalibrator - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar\
            -T PrintReads -nct 4 -R $HREFF -I patient3_n.sorted.dedup.bam -BQSR patient3_n.recal.table \
            -o patient_3_n.final.bam
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
            -T PrintReads -nct 4 -R $HREFF -I patient3_t.sorted.dedup.bam -BQSR patient3_t.recal.table \
            -o patient_3_t.final.bam


## Somatic mutation calling (BAM file -> VCF file)
Since we do not have time and capacity to run a whole sample during our exercises we will call somatic mutations on a single chromosome of your choice. Simply choose chromosome name before runnig mutect (e.g. CHR=chr15).
MuTect2 is a somatic mutation caller developed by Broad Institute. MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect [(Cibulskis et al., 2013)](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html) with the assembly-based machinery of [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php). The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller. To learn more about Mutect2 follow this link [MuTect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)

        ### Set chromosome:
        CHR=TYPE_CHROMOSOME_HERE
        ### Run Mutect2
        time java -Xmx4G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=1 
            -jar ${GATKROOT}/GenomeAnalysisTK.jar -T MuTect2 -R $HREFF --dbsnp $SREFF --cosmic $cosmicREFF \
            -I:tumor patient_3_t.final.bam -I:normal patient_3_n.final.bam -o patient_3_t.${CHR}.mutect2.vcf -L ${CHR}
        ### To run a whole genome simply do not use the -L option.

## Inference of tissue of origin

        ### .....

## Inference of copy number profile

        ### .....

## Data visualization

        ### .....

