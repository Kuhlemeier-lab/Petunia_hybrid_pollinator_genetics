#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="7combgvcf"

# Runtime and memory
#SBATCH --time=96:00:00  # 15h 
#SBATCH --mem-per-cpu=64G  # 14Gb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/07_combine_gvcfs_%A.out
#SBATCH --error=code/07_combine_gvcfs_%A.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

echo -e "#### Combine and then genotype gvcf files
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add Development/java/1.8.0_202
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0

# print which modules are loaded (gives the version of the software
# that is used
module list 2>&1
# print detail of java version used
java -showversion 2>&1 >/dev/null | grep -P "JDK|jdk"

# reference genome
genome=${scdir}/data/genomes/Peax403.fasta

# combine GVCFs
cd ${scdir}/data/raw/variants
echo -e "\n## Combine GVCFs \ncurrent dir: `pwd`\n"

GenomeAnalysisTK CombineGVCFs \
        -R ${genome} \
        --variant 1.g.vcf \
        --variant 2.g.vcf \
        --variant 3.g.vcf \
        --variant 4.g.vcf \
        --variant 5.g.vcf \
        --variant 6.g.vcf \
        --variant 7.g.vcf \
        --variant 8.g.vcf \
        --variant 9.g.vcf \
        --variant 10.g.vcf \
        --variant 11.g.vcf \
        --variant 12.g.vcf \
        --variant 13.g.vcf \
        --variant 14.g.vcf \
        --variant 15.g.vcf \
        --variant 16.g.vcf \
        --variant 17.g.vcf \
        --variant 18.g.vcf \
        --variant 19.g.vcf \
        --variant 20.g.vcf \
        --variant 21.g.vcf \
        --variant 22.g.vcf \
        --variant 23.g.vcf \
        --variant 24.g.vcf \
        --variant 25.g.vcf \
        --variant 26.g.vcf \
        --variant 27.g.vcf \
        --variant 28.g.vcf \
        --variant 29.g.vcf \
        --variant 30.g.vcf \
        --variant 31.g.vcf \
        --variant 32.g.vcf \
        --variant 33.g.vcf \
        --variant 34.g.vcf \
        --variant 35.g.vcf \
        --variant 36.g.vcf \
        --variant 37.g.vcf \
        --variant 38.g.vcf \
        --variant 39.g.vcf \
        --variant 40.g.vcf \
        --variant 41.g.vcf \
        --variant 42.g.vcf \
        --variant 43.g.vcf \
        --variant 44.g.vcf \
        --variant 45.g.vcf \
        --variant 46.g.vcf \
        --variant 47.g.vcf \
        --variant 48.g.vcf \
        --variant 49.g.vcf \
        --variant 50.g.vcf \
        --variant 51.g.vcf \
        --variant 52.g.vcf \
        --variant 53.g.vcf \
        --variant 54.g.vcf \
        --variant 55.g.vcf \
        --variant 56.g.vcf \
        --variant 57.g.vcf \
        --variant 58.g.vcf \
        --variant 59.g.vcf \
        --variant 60.g.vcf \
        --variant 61.g.vcf \
        --variant 62.g.vcf \
        --variant 63.g.vcf \
        --variant 64.g.vcf \
        --variant 65.g.vcf \
        --variant 66.g.vcf \
        --variant 67.g.vcf \
        --variant 68.g.vcf \
        --variant 69.g.vcf \
        --variant 70.g.vcf \
        -O cohort.g.vcf

echo -e "\n## Combine done."
