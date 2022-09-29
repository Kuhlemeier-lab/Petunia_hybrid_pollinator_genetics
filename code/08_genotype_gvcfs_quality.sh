#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="8geno"

# Runtime and memory
#SBATCH --time=48:00:00    # 15h 
#SBATCH --mem-per-cpu=10G  # 4 Gb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/08_genotype_gvcfs_quality_%A.out
#SBATCH --error=code/08_genotype_gvcfs_quality_%A.err
#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Genotype gvcf file
## `date`
## current dir: `pwd`
####\t####\n"

module load Workspace
module load vital-it
module load R/3.4.2
module add UHTS/Analysis/picard-tools/2.21.8
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module list 2>&1

# print which modules are loaded (gives the version of the software
# that is used)
module list 2>&1


# reference genome
genome=${scdir}/data/genomes/Peax403.fasta

cd ${scdir}/data/raw/variants

GenomeAnalysisTK GenotypeGVCFs \
     -R ${genome} \
     -V cohort.g.vcf \
     -O raw.vcf \
     --seconds-between-progress-updates 600

echo -e "\nGenotype done.\n Split into SNPs and INDELs for quality plotting.\n"

picard-tools SplitVcfs \
    I=raw.vcf \
    SNP_OUTPUT=raw_snp.vcf \
    INDEL_OUTPUT=raw_indel.vcf \
    STRICT=false

echo -e "\nExtracting a table for quality parameters plotting.\n"
GenomeAnalysisTK VariantsToTable \
     -R ${genome} \
     -V raw_snp.vcf \
     -O raw_snp.vcf.table \
     -F CHROM \
     -F POS \
     -F QD \
     -F FS \
     -F SOR \
     -F MQ \
     -F MQRankSum \
     -F ReadPosRankSum \
     --show-filtered

GenomeAnalysisTK  VariantsToTable \
    -R ${genome} \
    -V raw_indel.vcf \
    -O raw_indel.vcf.table \
    -F CHROM \
    -F POS \
    -F QD \
    -F FS \
    -F SOR \
    -F ReadPosRankSum \
    --show-filtered

echo -e "\nraw_snp.vcf.table and raw_indel.vcf.table were saved in `pwd`\n"

# call plot_distribution.R
echo -e "Make distribution plots of SNPs and INDELs.\n"

if [ ! -d ${chdir}/figures/exploratory ]; then
	mkdir -p ${chdir}/figures/exploratory
fi
# with GATK 4.1.3 I get a trailing tab on the header line of the table,
# so I remove it beforehand to avoid problems in R.
sed -i "s/\t$//" raw_indel.vcf.table
sed -i "s/\t$//" raw_snp.vcf.table
# requires optparse library
Rscript ${chdir}/code/plot_vcfq_distribution.R -i raw_indel.vcf.table --vcftype INDEL -o ${chdir}/figures/exploratory/
Rscript ${chdir}/code/plot_vcfq_distribution.R -i raw_snp.vcf.table --vcftype SNP -o ${chdir}/figures/exploratory/
echo -e "Plots done. They are saved in figures/exploratory/ .\n"
