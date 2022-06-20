#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="12pcangsd"

# Runtime and memory
#SBATCH --time=12:10:00  #1h
#SBATCH --mem-per-cpu=16G  #6Gb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/12a_pcangsd_%A.out
#SBATCH --error=code/12a_pcangsd_%A.err
#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

module load Anaconda3
# pcangsd version 1.10


echo -e "#### PCAngsd
## `date`
## current dir: `pwd`
####\t####\n"

eval "$(conda shell.bash hook)"
conda activate pcangsd

cd ${scdir}/data/raw/variants

if [ ! -d ${scdir}/data/raw/pca ]
then
    mkdir ${scdir}/data/raw/pca
fi

pcangsd -b ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.beagle.gz \
  -o ${scdir}/data/raw/pca/cr09_mm005 \
  --minMaf 0 \
  -t 1

