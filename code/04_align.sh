#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="4align"

# Runtime and memory
#SBATCH --time=24:00:00  # 2h
#SBATCH --mem-per-cpu=4G  # 20Gb
#SBATCH --cpus-per-task=16

#SBATCH --array=1-70

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/04_align_%A_%a.out
#SBATCH --error=code/04_align_%A_%a.err

#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

genome=${scdir}/data/genomes/Peax403.fasta

echo -e "#### Align reads
## `date`
## current dir: `pwd`
####\t####\n"

module add vital-it
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.10

module list 2>&1

# align reads to the P. axillaris unmasked 403 assembly
cd ${scdir}/data/raw
if [ ! -d aligned_reads ]; then
	mkdir aligned_reads
fi

cd ${scdir}/data/raw/rawreads

# The reads from the bioproject PRJNA522653 are already grouped in forward
# and reverse (no need to merge the lanes).
# The reads from PRJNA706535 are not merged, so I map them and then merge them.
# The samples that don't need merging are listed in the array.

# declare an array of indexes for the already merged samples
# needed for the if statement. Otherwise the list is too long.
declare -a mergedIDs=("2" "3" "5" "6" "8" "9"
                   "10" "11" "12" "20" "21"
                   "24" "25" "33" "34" "42"
                   "46" "48" "52" "53" "54" "57")

if [[ " ${mergedIDs[*]} " == *" ${SLURM_ARRAY_TASK_ID} "* ]];
then 
  bwa mem -t 16 \
    -M ${genome} \
    -R "@RG\tID:${SLURM_ARRAY_TASK_ID}\tSM:${SLURM_ARRAY_TASK_ID}\tLB:2016lib${SLURM_ARRAY_TASK_ID}\tPL:illumina\tDT:2016-09-01" \
    ${SLURM_ARRAY_TASK_ID}_L1_SRR*.sra_1.fastq.gz.trimmed.fastp.fq.gz ${SLURM_ARRAY_TASK_ID}_L1_SRR*.sra_2.fastq.gz.trimmed.fastp.fq.gz | samtools sort -@ 16 -o ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L1.bam -
else
  # lane 1
  bwa mem -t 16 -M ${genome} -R "@RG\tID:${SLURM_ARRAY_TASK_ID}L1\tSM:${SLURM_ARRAY_TASK_ID}\tLB:2018lib${SLURM_ARRAY_TASK_ID}\tPL:illumina\tDT:2018-09-01" ${SLURM_ARRAY_TASK_ID}_L1_SRR*.sra_1.fastq.gz.trimmed.fastp.fq.gz ${SLURM_ARRAY_TASK_ID}_L1_SRR*.sra_2.fastq.gz.trimmed.fastp.fq.gz | samtools sort -@ 16 -o ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L1.bam -
  
  # lane 2
  bwa mem -t 16 -M ${genome} -R "@RG\tID:${SLURM_ARRAY_TASK_ID}L2\tSM:${SLURM_ARRAY_TASK_ID}\tLB:2018lib${SLURM_ARRAY_TASK_ID}\tPL:illumina\tDT:2018-09-01" ${SLURM_ARRAY_TASK_ID}_L2_SRR*_1.fastq.gz.trimmed.fastp.fq.gz ${SLURM_ARRAY_TASK_ID}_L2_SRR*.sra_2.fastq.gz.trimmed.fastp.fq.gz | samtools sort -@ 16 -o ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L2.bam -

  # merge bam files from different lanes
  samtools merge ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}.bam ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L1.bam ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L2.bam

  rm ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L1.bam
  rm ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_L2.bam
fi 

echo -e "\n## Done aligning.\n"
