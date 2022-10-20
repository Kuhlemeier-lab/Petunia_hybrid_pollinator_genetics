#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="fst"

# Runtime and memory
#SBATCH --time=12:00:00  # 9h. 12' to do saf, 240Mb + 4h, 2gb with all.
#SBATCH --mem-per-cpu=8G  # 3Gb
#SBATCH --cpus-per-task=1

#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/14a_divergence_%A.out
#SBATCH --error=code/14a_divergence_%A.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

echo -e "#### Fst from ANGSD
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Analysis/samtools/1.10
angsddir=/xxx/tools/angsd      # version 0.933-111-g5859d2b

module list 2>&1

if [ ! -d ${scdir}/data/raw/fst ]
then
  mkdir -p ${scdir}/data/raw/fst
fi 

# reference genome
genome=${scdir}/data/genomes/Peax403.fasta

# using this dataset:
base_vcf=${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.recode.vcf
# and samples grouped as in file 
groups_df="${chdir}/data/raw/samples_in_admixture_groups.csv"

grouping_var="admixture"
echo -e "Using ${grouping_var} as grouping variable\n"
# and which group in that variable
group=1
echo -e "Group ${group}\n"
group_members=`grep -oP "^${grouping_var},${group},\K[0-9,]+" ${groups_df}`
echo -e "The group members are:\n${group_members}\n"
bcftools view -s "${group_members}"  ${base_vcf} > ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_${grouping_var}_${group}.vcf
# count samples
n_members1=`echo $group_members | awk '{split($1,a,",");for(i in a)print a[i]}' | grep -cP "^\d"`
echo -e "There is ${n_members1} members in the group\n"

# extract group 2 samples
group=2
group_members=`grep -oP "^${grouping_var},${group},\K[0-9,]+" ${groups_df}`
echo -e "Group ${group}\n" 
echo -e "The group members are:\n${group_members}\n" 
bcftools view -s "${group_members}" ${base_vcf} > ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_${grouping_var}_${group}.vcf
# count samples   
n_members2=`echo $group_members | awk '{split($1,a,",");for(i in a)print a[i]}' | grep -cP "^\d"` 
echo -e "There is ${n_members2} members in the group\n"  

# calculate saf (sample allele frequency) for group 1
${angsddir}/angsd \
  -vcf-PL ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_admixture_1.vcf \
  -doSaf 1 \
  -out ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_1 \
  -nInd ${n_memebers1} \
  -doMajorMinor 1 \
  -anc ${genome} 

# calculate saf for group 2
${angsddir}/angsd \
  -vcf-PL ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005_admixture_2.vcf \
  -doSaf 1 \
  -out ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2 \
  -nInd ${n_members2} \
  -doMajorMinor 1 \
  -anc ${genome}

# steps to check doSaf has worked properly
#realSFS print ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2.saf.idx > ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2.saf.idx.table
#realSFS check ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2.saf.idx 
#${angsdDir}/misc/realSFS check ${scdir}/data/raw/fst/

# calculate 2d sfs prior, unfolded
${angsddir}/misc/realSFS \
  -P 1 \
  -maxIter 50000 \
  -tole 0.00001 \
  ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_1.saf.idx \
  ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2.saf.idx > ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded.ml

# calculate fst, unfolded
${angsddir}/misc/realSFS fst index \
  ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_1.saf.idx \
  ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_2.saf.idx \
  -sfs ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded.ml \
  -fstout ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded \
  -P 1

# get global estimate
${angsddir}/misc/realSFS fst stats ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded.fst.idx -P 1
# 200 kb window, 100 kb step
${angsddir}/misc/realSFS fst stats2 ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded.fst.idx \
  -P 1 \
  -win 200000 \
  -step 100000 > ${scdir}/data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded_slidingwindow200kb100kb
