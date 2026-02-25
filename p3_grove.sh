#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cmg276@pitt.edu
#SBATCH -M teach
#SBATCH -A hugen2072-2026s
#SBATCH -t 48:00:00
#SBATCH --mem=64g
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH -N 1
#SBATCH -c 24

#load the modules needed
module load fastqc/0.11.9
module load cutadapt/5.1
module load gcc/8.2.0 bwa/0.7.17 samtools/1.22.1
module load gatk/4.5.0.0

#files are in /ix1/hugen2072-2026s/p3 - shortened to p3
#toy dataset = toy_{1,2}.fastq.gz

#Read Quality Control
#Skipped this step for p3


#Filter reads
fastqc p3/p3_1.fastq.gz -t 24 --outdir .
fastqc p3/p3_2.fastq.gz -t 24 --outdir .

#What is good and bad about the reads?

#Trim the reads
cutadapt -j 0 -m 10 -q 20 p3/p3_1.fastq.gz p3/p3_2.fastq.gz \
-a AGATCGGAAGAG -A AGATCGGAAGAG \
-o $SLURM_SCRATCH/p3_1_trimmed.fastq.gz -p $SLURM_SCRATCH/p3_2_trimmed.fastq.gz

#Check the trimmed Reads
fastqc $SLURM_SCRATCH/p3_1_trimmed.fastq.gz -t 24 --outdir=.
fastqc $SLURM_SCRATCH/p3_2_trimmed.fastq.gz -t 24 --outdir=.

#Alignment
bwa mem -t 24 \
p3/Homo_sapiens_assembly38.fasta \
$SLURM_SCRATCH/p3_1_trimmed.fastq.gz $SLURM_SCRATCH/p3_2_trimmed.fastq.gz \
-R "@RG\tID:p3\tLB:P5\tSM:P5\tPL:ILLUMINA" \
> $SLURM_SCRATCH/p3.sam
samtools view -b $SLURM_SCRATCH/p3.sam | samtools sort -o $SLURM_SCRATCH/p3.bam
samtools index $SLURM_SCRATCH/p3.bam

#put the bam files as Cram in home directory
samtools view -C -T p3/Homo_sapiens_assembly38.fasta \
  --output-fmt-option version=3.0 \
  -o p3.cram $SLURM_SCRATCH/p3.bam
samtools index p3.cram && \
rm $SLURM_SCRATCH/p3.sam $SLURM_SCRATCH/p3.bam



#Alignment Quality Control
gatk MarkDuplicatesSpark \
  -I p3.cram \
  -O $SLURM_SCRATCH/p3_dupsmarked.bam \
  -R p3/Homo_sapiens_assembly38.fasta

#Alignment Mismatches
#recalibration stats
gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=24" \
-I $SLURM_SCRATCH/p3_dupsmarked.bam \
-R p3/Homo_sapiens_assembly38.fasta \
-O $SLURM_SCRATCH/p3_BQSR.table \
--known-sites p3/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites p3/dbsnp_146.hg38.vcf.gz

#apply the recal stats to bam file
gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=24" \
-R p3/Homo_sapiens_assembly38.fasta \
-I $SLURM_SCRATCH/p3_dupsmarked.bam \
--bqsr-recal-file $SLURM_SCRATCH/p3_BQSR.table \
-O $SLURM_SCRATCH/p3_dupsmarked_cleaned.bam

#create another cram file
samtools view -C -T p3/Homo_sapiens_assembly38.fasta \
  --output-fmt-option version=3.0 \
  -o p3_dup_clean.cram $SLURM_SCRATCH/p3_dupsmarked_cleaned.bam
samtools index p3_dup_clean.cram

samtools flagstat -@ 24 p3.cram \
> p3_alignment_statistics_grove.out
samtools depth p3.cram \
| gzip > p3_depth_statistics.out.gz

#how many total reads passed QC in pipeline?


#how many were mapped and what % of QC-passed reads is that?


#what is the total read depth and what is the mapped read depth?






#Genotyping
#define active regions, determine haplotypes by reassembly of the active region,
#determine likilyhoods of the haplotypes given real data, and assign sample genotypes
gatk HaplotypeCaller -R p3/Homo_sapiens_assembly38.fasta \
-I p3.cram \
-O $SLURM_SCRATCH/p3_genotypes.g.vcf.gz \
-ERC GVCF \
-OVI \
--native-pair-hmm-threads 24


#convert thr gvcf to a vcf
gatk GenotypeGVCFs --java-options "-XX:ParallelGCThreads=24" \
-R p3/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/p3_genotypes.g.vcf.gz \
-O $SLURM_SCRATCH/p3_genotypes.vcf.gz





#Genotype Quality Control

#vcf with site info
gatk MakeSitesOnlyVcf -I $SLURM_SCRATCH/p3_genotypes.vcf.gz -O $SLURM_SCRATCH/p3_sites_only.vcf.gz


#find INDEL sites and recalibrate the genotype
gatk VariantRecalibrator --java-options "-XX:ParallelGCThreads=24" \
-mode INDEL \
-R p3/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/p3_sites_only.vcf.gz \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
-resource:mills,known=false,training=true,truth=true,prior=12 \
p3/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2 \
p3/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SCRATCH/p3_indels.recal \
--tranches-file $SLURM_SCRATCH/p3_indels.tranches


#table of recalibration stats
gatk VariantRecalibrator --java-options "-XX:ParallelGCThreads=24" \
-mode SNP \
-R p3/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/p3_sites_only.vcf.gz \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-resource:hapmap,known=false,training=true,truth=true,prior=15 \
p3/hapmap_3.3.hg38.vcf.gz \
-resource:omni,known=false,training=true,truth=true,prior=12 \
p3/1000G_omni2.5.hg38.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10 \
p3/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=7 \
p3/dbsnp_146.hg38.vcf.gz \
-O $SLURM_SCRATCH/p3_snps.recal \
--tranches-file $SLURM_SCRATCH/p3_snps.tranches


#gatk ApplyVQSR --java-options
#used in steps 12c and 12 d instead
#apply recalibration stats
gatk ApplyVQSR --java-options "-XX:ParallelGCThreads=24" \
-mode INDEL \
-R p3/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/p3_genotypes.vcf.gz \
--recal-file $SLURM_SCRATCH/p3_indels.recal \
--tranches-file $SLURM_SCRATCH/p3_indels.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-O $SLURM_SCRATCH/p3_genotypes_indelqc.vcf.gz

#apply recal stats to SNPs
gatk ApplyVQSR --java-options "-XX:ParallelGCThreads=24" \
-mode SNP \
-R p3/Homo_sapiens_assembly38.fasta \
-V $SLURM_SCRATCH/p3_genotypes_indelqc.vcf.gz \
--recal-file $SLURM_SCRATCH/p3_snps.recal \
--tranches-file $SLURM_SCRATCH/p3_snps.tranches \
--truth-sensitivity-filter-level 99.0 \
--create-output-variant-index true \
-O p3_grove.vcf.gz


#Calculate some stats on the genotypes calls to assess quality
gatk CollectVariantCallingMetrics -I p3_grove.vcf.gz \
--DBSNP p3/dbsnp_146.hg38.vcf.gz \
-O p3_genotype_metrics_grove


#How many SNPs were called?


#How many indels were called?


#What is the Ts/Tv (transition/transversion) ratio for this call set?


#What is the indel (insertion/deletion) ratio for this call set?


#Are these statistics consistent with what we would expect for sequencing?


#Do you think these reads are from whole-exome or whole-genome sequencing?

























































