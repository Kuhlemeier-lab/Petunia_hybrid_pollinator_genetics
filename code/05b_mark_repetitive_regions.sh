#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="5repet"

# Runtime and memory
#SBATCH --time=48:10:00  # 3d 
#SBATCH --mem-per-cpu=16G  # 12Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/05b_mark_repetitive_regions_%A.out
#SBATCH --error=code/05b_mark_repetitive_regions_%A.err

#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Mark repetitive regions in the genome with repeat modeler
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module load SequenceAnalysis/Repeat/RepeatModeler/1.0.11
module list 2>&1

genome=${scdir}/data/genomes/Peax403.fasta

if [ ! -d ${scdir}/data/genomes/repeatModeler ]; then
    mkdir -p ${scdir}/data/genomes/repeatModeler
fi

cd ${scdir}/data/genomes/repeatModeler


# build a database
BuildDatabase \
    -name Peax403_repmodeler \
    -engine ncbi \
    ${genome}

# model repeats
RepeatModeler \
    -database Peax403_repmodeler \
    -engine ncbi \
    -pa 8

#downloaded TE database TREP from UZH accessed on 12/04/2022
# release 19, complete TREP nucleotide sequences
# This database was combined with the P. axillaris repeats
# just modelled above.

# mask repeats
RepeatMasker \
    -pa 4 \
    -s \
    -e ncbi \
    -lib ${scdir}/data/genomes/repeatModeler/TREP_and_Peax403_repeat_library.fasta \
    -dir ${scdir}/data/genomes/repeatModeler \
    -gff \
    ${genome}
# rename repeatmasker output to a meaningful name
mv ${scdir}/data/genomes/repeatModeler/${genome}.out.gff ${scdir}/data/genomes/repeatModeler/${genome}_repeatModeler_TREP_repeats.gff
