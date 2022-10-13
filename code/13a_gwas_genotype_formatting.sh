#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="formatgemma"

# Runtime and memory
#SBATCH --time=12:00:00  # 1h 
#SBATCH --mem-per-cpu=8G  # 4Gb
#SBATCH --cpus-per-task=1

#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/13a_gwas_genotype_formatting_%A.out
#SBATCH --error=code/13a_gwas_genotype_formatting_%A.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids


module load vital-it
module add UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/vcftools/0.1.15

module list 2>&1

echo -e "#### Convert vcf GL to bimbam format
## `date`
## current dir: `pwd`
####\t####\n"

echo -e "filter from hardfiltered biallelic to cr09 mm005 but keep info field"
vcftools --gzvcf ${scdir}/data/raw/variants/hardfiltered_biallelic.vcf --positions ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.pos --recode --recode-INFO-all --out ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_vcfinfo

echo -e "hardfiltered_biallelic_cr09_mm005.pos has this many positions:\n"
wc -l ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.pos

echo -e "hardfiltered_biallelic_cr09_mm005_vcfinfo.recode.vcf has this many positions:\n"
grep -cPv "^#" ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_vcfinfo.recode.vcf

echo -e "Transform the vcf into a bimbam format file\n"
${chdir}/code/bcf2bbgeno_edit.pl -i ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_vcfinfo.recode.vcf -p H-W -o ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.bimbam -s -r
echo -e "The ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.bimbam file has this many positions:\n"
wc -l ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.bimbam

gzip ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.bimbam

