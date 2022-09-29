#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="6vars"

# Runtime and memory
#SBATCH --time=96:00:00  # 17 h with 32cpus
#SBATCH --mem-per-cpu=2G  # 13 Gb on 32cpus
#SBATCH --cpus-per-task=32

# Array
##SBATCH --array=2-70%20

# running tasks that require longer runtime
#SBATCH --array=46

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/06_call_variants_%A_%a.out
#SBATCH --error=code/06_call_variants_%A_%a.err

#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Call variants
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add Development/java/1.8.0_202
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0

# print which modules are loaded (gives the version of the software
# that is used)
module list 2>&1
# print detail of java version used
java -showversion 2>&1 >/dev/null | grep -P "JDK|jdk"

# reference genome
genome=${scdir}/data/genomes/Peax403.fasta

cd ${scdir}/data/raw/
if [ ! -d variants ]; then
  mkdir variants
fi
echo -e "\n## Call variants \ncurrent dir: `pwd`\n"

GenomeAnalysisTK HaplotypeCaller \
     -R ${genome} \
     -I ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam \
     -ERC GVCF \
     --do-not-run-physical-phasing \
     -O ${scdir}/data/raw/variants/${SLURM_ARRAY_TASK_ID}.g.vcf \
     --native-pair-hmm-threads 32 \
     --exclude-intervals ${scdir}/data/genomes/Peax403_repetitive_and_min100Xcov.bed

echo -e "Done variants calling.\n"
