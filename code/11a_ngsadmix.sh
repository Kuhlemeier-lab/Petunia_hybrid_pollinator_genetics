#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="11ngsadm"

# Runtime and memory
#SBATCH --time=96:00:00  #13h
#SBATCH --mem-per-cpu=48G  # 18Gb
#SBATCH --cpus-per-task=2 
#SBATCH --array=1-10  

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/11a_ngsadmix_%A_%a.out
#SBATCH --error=code/11a_ngsadmix_%A_%a.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

echo -e "#### NGSadmix 
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
angsddir=/xxx/tools/angsd
module list 2>&1

if [ ! -d ${scdir}/data/raw/admixture ]
then
  mkdir -p ${scdir}/data/raw/admixture
fi 


# Run for K {1-8}

for k in {1..8}
do
    ##################
    # I set the seed myself. It's based on the date in nanoseconds and
    # on the iteration number and K number.
    ##################
    # take date in nanosecods
    mydate=`date +%N`
    # set custom value
    myvalue=$k
    # add myvalue to mydate
    seed1=`echo "${mydate} + ${myvalue}" | bc`
    # concatenate myvalue after seed1 to pad in case a long
    # seed number is required
    seed2=${seed1}${myvalue}${myvalue}${myvalue}${myvalue}${myvalue}
    # crop seed at desired length (here 10 characters)
    seed3=${seed2:0:10}
    ##################
    # -P number of threads
    ${angsddir}/misc/NGSadmix -likes ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.beagle.gz -minMaf 0 -seed ${seed3} -K ${k} -P 2 -o ${scdir}/data/raw/admixture/hf_ba_cr09_mm005_k${k}_run${SLURM_ARRAY_TASK_ID}
done
