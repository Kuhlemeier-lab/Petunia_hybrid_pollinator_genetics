#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="5md"

# Runtime and memory
#SBATCH --time=24:10:00  # 3h  coverage takes about 1 h
#SBATCH --mem-per-cpu=12G  # 9Gb
#SBATCH --cpus-per-task=4

# Array
#SBATCH --array=1-70

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/05a_markdup_metrics_coverage_%A_%a.out
#SBATCH --error=code/05a_markdup_metrics_coverage_%A_%a.err

#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Mark duplicates and index bam files, calculate some stats, get regions with coverage > 100
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module add UHTS/Analysis/picard-tools/2.21.8
module add UHTS/Analysis/BEDTools/2.29.2

module list 2>&1

genome=${scdir}/data/genomes/Peax402INV.fasta

# Mark the duplicated reads
cd ${scdir}/data/raw/aligned_reads

if [ ! -d temp_gatk ]; then
    mkdir temp_gatk
fi

if [ -f ${SLURM_ARRAY_TASK_ID}_L1.bam ]; then
  bamfile=${SLURM_ARRAY_TASK_ID}_L1.bam
else
  bamfile=${SLURM_ARRAY_TASK_ID}.bam
fi

GenomeAnalysisTK MarkDuplicatesSpark \
    -I ${bamfile} \
    -O ${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam \
    -M ${SLURM_ARRAY_TASK_ID}_marked_dup_metrics.txt \
    --duplicate-tagging-policy "All" \
    --create-output-bam-index "true" \
    --conf 'spark.executor.cores=4' \
    --conf "spark.local.dir=temp_gatk" \
    --verbosity WARNING

# remove old bam file
rm ${SLURM_ARRAY_TASK_ID}.bam

echo -e "Done duplicates. Calculate more metrics...\n\n"

picard-tools CollectWgsMetrics \
    I=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam \
    O=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates_wgsmetrics.txt \
    R=${genome}

echo -e "calculate wgs metris only on gene space"
# I have to make a picard interval list file
# first I get the header for the file
grep -P "(@SQ)|(@HD)" ${scdir}/data/genomes/Peax402INV.dict > ${scdir}/data/genomes/Peax402INV.gene${SLURM_ARRAY_TASK_ID}.picardlist
# extract gene regions into a picard interval list file

ln -s /xxx/peaxi162AQ_Peax402INV.cds.gff ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff

grep -P "\sgene\s" ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff | cut -f 1,4,5,7,9 | cut -d ";" -f 1 >> ${scdir}/data/genomes/Peax402INV.gene${SLURM_ARRAY_TASK_ID}.picardlist
# calculate metrics on this subset of intervals
picard-tools CollectWgsMetrics \
    I=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam \
    O=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates_genespace_wgsmetrics.txt \
    R=${genome} \
    INTERVALS=${scdir}/data/genomes/Peax402INV.gene${SLURM_ARRAY_TASK_ID}.picardlist

rm ${scdir}/data/genomes/Peax402INV.genesonly.gene${SLURM_ARRAY_TASK_ID}.picardlist

picard-tools CollectQualityYieldMetrics \
    I=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam \
    O=${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates_qualityyieldmetrics.txt


export SINGULARITY_BINDPATH="$scdir:/scdir"
singularity exec /storage/homefs/mbinaghi/covtobed_latest.sif coverage -q -s -m 100 -r -i ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_marked_duplicates.bam > ${scdir}/data/raw/aligned_reads/${SLURM_ARRAY_TASK_ID}_md_mincov100
