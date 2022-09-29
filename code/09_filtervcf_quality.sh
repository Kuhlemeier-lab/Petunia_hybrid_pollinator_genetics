#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="9filt"

# Runtime and memory
#SBATCH --time=48:00:00  # 2h 
#SBATCH --mem-per-cpu=24G  # 400Mb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/hybrids_peaxiINV
#SBATCH --output=code/09_filtervcf_quality_%A.out
#SBATCH --error=code/09_filtervcf_quality_%A.err
#################################
chdir=/xxx/hybrids_peaxiINV
scdir=/xxx/hybrids_peaxiINV

echo -e "#### Filter and select relevant variants.

## Apply hard filters to SNPs and INDELs separately.
##  `date`
##  current dir: `pwd`
####\t####\n
"
module load vital-it
module add UHTS/Analysis/picard-tools/2.21.8
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/vcftools/0.1.15
angsddir=/xxx/tools/angsd


# print which modules are loaded (gives the version of the software
# that is used)
module list 2>&1
# print detail of java version used
java -showversion 2>&1 >/dev/null | grep -P "JDK|jdk"

# reference genome
genome=${scdir}/data/genomes/Peax403.fasta

# Filter variants
cd ${scdir}/data/raw/variants

echo -e "How many variants in the raw set? "
grep -Pcv "^#" raw.vcf

echo -e "## Hard filter SNP variants
"

GenomeAnalysisTK VariantFiltration \
  --add-output-vcf-command-line \
  -R ${genome} \
  -V raw_snp.vcf \
  --filter-expression "QD<2.0" \
  --filter-name "QD" \
  --filter-expression "FS>60.0" \
  --filter-name "FS" \
  --filter-expression "SOR>3.0" \
  --filter-name "SOR" \
  --filter-expression "MQ<40.0" \
  --filter-name "MQ" \
  --filter-expression "MQRankSum<-12.5" \
  --filter-name "MQRankSum" \
  --filter-expression "ReadPosRankSum<-8.0" \
  --filter-name "ReadPosRankSum" \
  -O hardfiltered_snp.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
  --add-output-vcf-command-line \
  -R ${genome} \
  -V hardfiltered_snp.vcf \
  --exclude-filtered \
  -O hf_selected_snp.vcf

echo -e "Select done.
## Hard filter INDEL variants
"

GenomeAnalysisTK VariantFiltration \
  --add-output-vcf-command-line \
  -R ${genome} \
  -V raw_indel.vcf \
  --filter-expression "QD<2.0" \
  --filter-name "QD" \
  --filter-expression "FS>200.0" \
  --filter-name "FS" \
  --filter-expression "SOR>10.0" \
  --filter-name "SOR" \
  --filter-expression "ReadPosRankSum<-20.0" \
  --filter-name "ReadPosRankSum" \
  -O hardfiltered_indel.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
  --add-output-vcf-command-line \
  -R ${genome} \
  -V hardfiltered_indel.vcf \
  --exclude-filtered \
  -O hf_selected_indel.vcf

# merge SNP and INDEL variants in a single file
echo -e "Selection done.
Merging SNPs and INDELs...
"

picard-tools MergeVcfs \
  I=hf_selected_snp.vcf \
  I=hf_selected_indel.vcf \
  O=hardfiltered.vcf

echo -e "\n How many variants in the hardfiltered.vcf set? "
grep -Pcv "^#" hardfiltered.vcf


echo -e "Hard filter, select, merge SNP and INDEL done.
Selecting biallelic.
"

GenomeAnalysisTK SelectVariants \
  --add-output-vcf-command-line \
  -R ${genome} \
  -V hardfiltered.vcf \
  --restrict-alleles-to BIALLELIC \
  -O hardfiltered_biallelic.vcf

echo -e "How many variants in the hardfiltered_biallelic.vcf set? "
grep -Pcv "^#" hardfiltered_biallelic.vcf


echo -e "#### Make a table with AF and DP values for plotting"

GenomeAnalysisTK VariantsToTable \
  -R ${genome} \
  -V hardfiltered_biallelic.vcf \
  -O hardfiltered_biallelic_AFDPMQ.vcf.table \
  -F CHROM \
  -F POS \
  -F AF \
  -F DP \
  -F MQ

echo -e "#### Filter variants for call rate and MAF, select only SNPs"

vcftools --vcf ${scdir}/data/raw/variants/hardfiltered_biallelic.vcf \
  --remove-indels \
  --max-missing 0.9 \
  --remove-filtered-geno-all \
  --recode \
  --out ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09

echo -e "Count positions with call rate 0.90 in hardfiltered biallelic SNPs: \n"
grep -vc "^#" ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09.recode.vcf

echo -e "Use ANGSD to get a list of positions with MAF > 0.05\n"
# note that some versions of ANGSD tools mess up the output file when 
# used in multithreaded
${angsddir}/angsd -doMaf 1 \
  -doMajorMinor 1 \
  -out ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005 \
  -vcf-PL ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09.recode.vcf \
  -nInd 70 -P 1 \
  -fai ${scdir}/data/genomes/Peax403.fasta.fai \
  -minMaf 0.05

echo -e "Re-format mafs file to use it to subset the vcf with vcftools\n"
gunzip -c ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.mafs.gz | grep -Pv "^chromo" | cut -f 1,2 > ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.pos

echo -e "This many SNPs are in the hardfiltered_biallelic_cr09_mm005.pos  "
wc -l hardfiltered_biallelic_cr09_mm005.pos

echo -e "Subset the vcf to keep only positions with maf > 0.05\n"
vcftools --vcf ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09.recode.vcf \
  --positions ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005.pos \
  --recode \
  --out ${scdir}/data/raw/variants/hardfiltered_biallelic_cr09_mm005

echo -e "## Filters done.\n"
