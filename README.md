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

### 1.1 Take a first look at the data




## 2. Data Pre-processing

### 2.1 Cleanup and alignment (FASTQ file -> BAM file)

1. Step 1 - (PLEASE DO NOT RUN AND GO DIRECTLY TO STEP 2)
Use bwa mem for aligning reads. Tumor sample and normal separately.
 
        ### Define variables for convenience
        ReadGoupID_N="@RG\tID:TCRBOA3-N-WEX\tSM:TCRBOA3\tPL:ILLUMINA\tLB:libN\tPU:TCRBOA3-N-WEX"
        ReadGoupID_T="@RG\tID:TCRBOA3-T-WEX\tSM:TCRBOA3\tPL:ILLUMINA\tLB:libT\tPU:TCRBOA3-T-WEX"
        f1n="/home/ngscourse/stud059/test_server1/somatic_calling/TCRBOA3-N-WEX.read1.fastq.bz2"
        f2n="/home/ngscourse/stud059/test_server1/somatic_calling/TCRBOA3-N-WEX.read2.fastq.bz2"
        f1t=TCRBOA3-T-WEX.read1.fastq.bz2
        f2t=TCRBOA3-T-WEX.read2.fastq.bz2
        HREFF="/home/ngscourse/stud059/test_server1/somatic_calling/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set"
        IREFF="/home/ngscourse/stud059/test_server1/somatic_calling/human_GRCh38/Indel_refs/mills_gold.b38.vcf"
        SREFF="/home/ngscourse/stud059/test_server1/somatic_calling/human_GRCh38/SNP_refs/1000G.snps.b38.vcf"
        cosmicREFF="/home/ngscourse/stud059/test_server1/somatic_calling/human_GRCh38/cosmic/CosmicCodingMuts_chr_sorted.vcf"
        GATKROOT=/home/ngscourse/stud059/test_server1/snp_calling
        
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

        PICARD=/home/ngscourse/stud059/test_server1/somatic_calling/picard-2.jar
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
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_n.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient3_n.recal.table
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
            -T BaseRecalibrator -nct 4 -R $HREFF -I patient3_t.sorted.dedup.bam -knownSites $SREFF \
            -knownSites $IREFF -o patient3_t.recal.table
      
        
6. Step 6 - BaseRecalibrator - Part 2.

        ### Run Base Recalibrator 2 - parallelize by using -nct option
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
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
            -o $patient_3_n.final.bam
        java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=4 -jar ${GATKROOT}/GenomeAnalysisTK.jar \
            -T PrintReads -nct 4 -R $HREFF -I patient3_t.sorted.dedup.bam -BQSR patient3_t.recal.table \
            -o $patient_3_t.final.bam


## Somatic mutation calling (BAM file -> VCF file)
Since we do not have time and capacity to run a whole sample during our exercises we will call somatic mutations on a single chromosome of your choice. Simply choose chromosome name before runnig mutect (e.g. CHR=chr15)

        ### Set chromosome:
        CHR=TYPE_CHROMOSOME_HERE
        ### Run Mutect2
        time java -Xmx10G -Xms1024M -XX:+UseParallelGC -XX:ParallelGCThreads=1 
            -jar ${GATKROOT}/GenomeAnalysisTK.jar -T MuTect2 -R $HREFF --dbsnp $SREFF --cosmic $cosmicREFF \
            -I:tumor $tumor -I:normal $normal -o $out_dir/${fnam}.chr${CHR}.mutect2.vcf -L chr${CHR}

## Inference of copy number profile

## Data visualization

