#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="gemma"

# Runtime and memory
#SBATCH --time=24:00:00  # 24h
#SBATCH --mem-per-cpu=48G  # 27Gb
#SBATCH --cpus-per-task=1

#SBATCH --partition=epyc2

#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids
#SBATCH --output=code/13c_gwas_bslmm_%A.out
#SBATCH --error=code/13c_gwas_bslmm_%A.err
#################################
chdir=/xxx/hybrids
scdir=/xxx/hybrids

# use this local version of gemma:
gemma=/xxx/tools/gemma-0.98.4-linux-static-AMD64

module list 2>&1

echo -e "#### Run association analysis
## `date`
## current dir: `pwd`
####\t####\n"

if [ ! -d ${scdir}/data/raw/gwas ]
then
  mkdir -p ${scdir}/data/raw/gwas
fi 

cd ${scdir}/data/raw/gwas/

dataset=hardfiltered_biallelic_cr09_mm005
# phenotype file
pheno_file=${scdir}/data/raw/gwas/pheno_gwas.bimbam
gt_data=${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.bimbam.gz

#------------------ Kinship matrix --------------------#
echo -e "calculate kinship matrix \n"
$gemma -g ${gt_data} \
       -p ${pheno_file} -gk 1 \
       -o ${dataset}_gwas \
       -hwe 0.001 -miss 0 -maf 0

#---------------------GWAS, LM -----------------------#
# lm p value method (1 = Wald, 2 = likelihood ratio test, 3 = score test, 4 = all)
pval=4
echo -e "\nLM\n"
# LM with no covariates:
for i in {1..3}; do
    pheno_idx=$i
    echo "phenotype ${pheno_idx}"
    echo ${dataset}
    gt_data=${gt_data}
    $gemma -g ${gt_data} \
           -p ${pheno_file} \
           -n ${pheno_idx} \
           -o ${dataset}_pheno${pheno_idx}_lm${pval} \
           -lm ${pval} \
           -hwe 0.001 \
           -miss 0 \
           -maf 0
done


#---------------------GWAS, LMM-----------------------#
# lmm p value method (1 = Wald, 2 = likelihood ratio test, 3 = score test, 4 = all)
pval=4
echo -e "\nLMM\n"
# LMM with no covariates:
for i in {1..3}; do
    pheno_idx=$i
    echo "phenotype ${pheno_idx}"
    echo ${dataset}
    gt_data=${gt_data}
    kinship=${scdir}/data/raw/gwas/output/${dataset}_gwas.cXX.txt
    $gemma -g ${gt_data} \
           -p ${pheno_file} \
           -n ${pheno_idx} \
           -o ${dataset}_pheno${pheno_idx}_lmm${pval} \
           -lmm ${pval} \
           -k ${kinship} \
           -hwe 0.001 \
           -miss 0 \
           -maf 0
done

#---------------------GWAS, BSLMM-----------------------#
echo -e "\nBSLMM\n"
dataset="hardfiltered_biallelic_cr09_mm005"
gt_data=${scdir}/data/raw/variants/${dataset}.bimbam.gz
pheno_file=${scdir}/data/raw/gwas/pheno_gwas.bimbam

# standard linear BSLMM 
# loop through phenotypes
for i in {1..3}; do
    pheno_idx=${i}
    echo "phenotype ${pheno_idx}"
    model_idx=1
    echo "Model -bslmm ${model_idx}"
    $gemma -g ${gt_data} \
           -p ${pheno_file} \
           -n ${pheno_idx} \
           -o ${dataset}_pheno${pheno_idx}_bslmm${model_idx} \
           -bslmm ${model_idx} \
           -hwe 0.001 \
           -miss 0 \
           -maf 0 \
           -w 50000000 \
           -s 200000000 \
           -rpace 1000 \
           -wpace 10000000
done

