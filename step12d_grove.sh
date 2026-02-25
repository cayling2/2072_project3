#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cmg276@pitt.edu
#SBATCH -M teach
#SBATCH -A hugen2072-2026s

module load gcc/15.1.0 bcftools/1.22 vcftools/0.1.16

#get the nonheader rows of the output
#use the -H tag found in slides which drops the header from the output

bcftools view p3_grove.vcf.gz -H | head -n 10 > grove_12d.txt



