#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="3index"

# Runtime and memory
#SBATCH --time=10:00:00  # 
#SBATCH --mem-per-cpu=32G  # 
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/03_index_dictionary_%A.out
#SBATCH --error=code/03_index_dictionary_%A.err

#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Make genome index
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/picard-tools/2.21.8
module add UHTS/Analysis/samtools/1.10

module list 2>&1

# make genome indexes for future use
ln -s /yyy/genomes/Peax402INV.fasta ${scdir}/data/genomes/Peax402INV.fasta
cd ${scdir}/data/genomes/
bwa index Peax402INV.fasta
echo -e "Indexing done.\n"

echo -e "Make dictionary.\n"
# Make picard genome dictionary for subsequent use
cd ${scdir}/data/genomes
echo -e "\n## dictionary\ncurrent dir: `pwd`\n"
picard-tools CreateSequenceDictionary \
    R=Peax402INV.fasta \
    O=Peax402INV.dict

samtools faidx Peax402INV.fasta
echo -e "Dictionary done.\n"

