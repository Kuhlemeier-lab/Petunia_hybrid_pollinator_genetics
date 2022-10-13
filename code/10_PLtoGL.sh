#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="pltogl"

# Runtime and memory
#SBATCH --time=12:00:00  # 2h 
#SBATCH --mem-per-cpu=48G  # 30Gb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/10_PLtoGL_%A.out
#SBATCH --error=code/10_PLtoGL_%A.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

echo -e "#### Convert PL values to Beagle GL
##  `date`
##  current dir: `pwd`
####\t####\n
"
module load vital-it

echo -e "Re-format hardfiltered_biallelic_cr09_mm005.recode.vcf\n"
python3 ${chdir}/code/gatkPLtobeagleGL.py -i ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.recode.vcf -o ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.beagle

sed -i "s/,/\\t/g" ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.beagle

# I run the file through dos2unix to remove weird characters that gave me troubles
# with ngsLD
gzip ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.beagle

