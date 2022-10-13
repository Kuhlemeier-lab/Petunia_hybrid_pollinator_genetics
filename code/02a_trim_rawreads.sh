#!/bin/bash
 
#SBATCH --mail-user=xxx
 
# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL
 
#SBATCH --job-name="2_trim"
 
# Runtime and memory
#SBATCH --time=06:00:00  # 35'
#SBATCH --mem-per-cpu=12G  # 4Gb
#SBATCH --cpus-per-task=12
 
# Array
#SBATCH --array=1-70

#SBATCH --partition=epyc2

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/02a_trim_rawreads_%A_%a.out
#SBATCH --error=code/02a_trim_rawreads_%A_%a.err
 
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids
 
echo -e "#### trimming
## `date`
## current dir: `pwd`
####\t####\n"
 
module load vital-it
module load UHTS/Quality_control/fastqc/0.11.7
module load UHTS/Analysis/trimmomatic/0.36
module load UHTS/Quality_control/fastp/0.19.5

module list 2>&1
 
# assess quality of each raw read file with fastqc
# 1 is the forward reads, 2 the reverse
echo -e "\n\n## Assessing fastq quality..."
cd ${scdir}/data/raw/rawreads
 
for file in ${SLURM_ARRAY_TASK_ID}_L*_SRR*.fastq.gz
do
    fastqc -t 12 ${file}
done
 
echo -e "\n\n## Run trimmomatic on lane 1..."
# for each file pair in the folder, run trimmomatic
for f1 in ${SLURM_ARRAY_TASK_ID}_L1_SRR*_1.fastq.gz
do
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
    trimmomatic PE -threads 12 -phred33 $f1 $f2 "$f1".trimmed.gz "$f1".unpaired.gz "$f2".trimmed.gz "$f2".unpaired.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

# remove polyG reads
    fastp -i "$f1".trimmed.gz -I "$f2".trimmed.gz -o "$f1".trimmed.fastp.fq.gz -O "$f2".trimmed.fastp.fq.gz -Q -L --trim_poly_g -h ${SLURM_ARRAY_TASK_ID}_L1 -w 12

done

echo -e "\n\n## Run trimmomatic on lane 2 if it exists"
for f1 in ${SLURM_ARRAY_TASK_ID}_L2_SRR*_1.fastq.gz
do
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
    trimmomatic PE -threads 12 -phred33 $f1 $f2 "$f1".trimmed.gz "$f1".unpaired.gz "$f2".trimmed.gz "$f2".unpaired.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

# remove polyG reads
    fastp -i "$f1".trimmed.gz -I "$f2".trimmed.gz -o "$f1".trimmed.fastp.fq.gz -O "$f2".trimmed.fastp.fq.gz -Q -L --trim_poly_g -h ${SLURM_ARRAY_TASK_ID}_L1 -w 12

done

echo -e "\n\n## Assess cleaned fastq quality"
for file in ${SLURM_ARRAY_TASK_ID}_*.fastq.gz.trimmed.fastp.fq.gz
do
    fastqc -t 12 ${file}
done

