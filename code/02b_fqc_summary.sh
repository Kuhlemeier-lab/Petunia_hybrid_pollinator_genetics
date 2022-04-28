#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="fqsum"

# Runtime and memory
#SBATCH --time=00:10:00  # 1'
#SBATCH --mem-per-cpu=4G  # 7Mb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/02b_fqc_summary_%A.out
#SBATCH --error=code/02b_fqc_summary_%A.err

#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Summary of fastqc outputs
## `date`
## current dir: `pwd`
####\t####\n"

module list 2>&1

python3 ${chdir}/code/fastqc_parser.py -p ${scdir}/data/raw/rawreads
