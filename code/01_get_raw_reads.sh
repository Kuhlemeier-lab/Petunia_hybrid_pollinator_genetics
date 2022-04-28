#!/bin/bash
 
#SBATCH --mail-user=xxx
 
# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL
 
#SBATCH --job-name="1_getsra"
 
# Runtime and memory
#SBATCH --time=00:10:00  # very long, especially the gzip step. Took 24h for about 20 entries.
#SBATCH --mem-per-cpu=4G  # very little, few Mb 
#SBATCH --cpus-per-task=6
 
# Partition
#SBATCH --partition=epyc2
#SBATCH --qos=job_epyc2_debug

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/01_get_raw_reads_%A.out
#SBATCH --error=code/01_get_raw_reads_%A.err
 
#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### get raw reads from SRA
## `date`
## current dir: `pwd`
####\t####\n"
 
module load vital-it
module add UHTS/Analysis/sratoolkit/2.10.7
module list 2>&1

cd ${chdir}

if [ ! -d ${scdir}/data/raw/rawreads ]; then
    mkdir -p ${scdir}/data/raw/rawreads
fi

cut -f 3 -d ";" ${chdir}/data/clean/reads_sample_ID.csv | grep -Po "SRR\d+" > ${scdir}/data/raw/rawreads/sra_ids.tmp

prefetch $(<${scdir}/data/raw/rawreads/sra_ids.tmp) --output-directory ${scdir}/data/raw/rawreads

fasterq-dump -e 6 --split-files --outdir ${scdir}/data/raw/rawreads ${scdir}/data/raw/rawreads/SRR*/*.sra

gzip ${scdir}/data/raw/rawreads/*.fastq

rm ${scdir}/data/raw/rawreads/SRR*/SRR*.sra
rm ${scdir}/data/raw/rawreads/sra_ids.tmp

# rename the files based on a csv
code/rename_reads.sh ${chdir}/data/clean/reads_sample_ID.csv ${scdir}/data/raw/rawreads

