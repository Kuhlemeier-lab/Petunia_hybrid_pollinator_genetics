# Petunia axillaris x P. exserta natural hybrids, genetic basis of a pollination syndrome

Bioinformatic analyses performed in the experiments for the publication about the genetic basis of pollination syndromes in Petunia axillaris and P. exserta.

Project conducted in the [Plant genetics and development group](https://www.ips.unibe.ch/research/deve/index_eng.html) at the Institute of Plant Sciences of the University of Bern.

Authors: Marta Binaghi, Marius Roesti, Therese Mandel, Korinna Esfeld, Loreta B. Freitas, Cris Kuhlemeier

Author of this page: Marta Binaghi

## DNA sequencing

Library preparation and sequencing were performed by the [Next Generation Sequencing platform of the University of Bern](https://www.ngs.unibe.ch/) in two batches, one in 2016 and one in 2018. In both batches, the DNA was amplified with illustraTM GenomiPhiTM V2 DNA Amplification Kit. Library preparation followed the TruSeq DNA PCR-free protocol. Sequencing was performed to obtain 150 bp long, paired-end reads, for an estimated coverage of 4-5 ×, on a genome size of 1.2 Gb (Bombarely et al., 2016). The 2016 batch was sequenced on two lanes of an Illumina HiSeq 3000. The 2018 batch was sequenced on two lanes (one chip) of the Illumina NovaSeq.

Raw reads were uploaded to NCBI SRA under [BioProjects PRJNA522653 (2016 batch)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA522653) and [PRJNA706535 (2018 batch)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA706535).

## Phenotype data

The phenotype of the plants was collected and is available in file [phenotype_sequenced_individuals.csv](data/phenotype_sequenced_individuals.csv). For more details on how the phenotypic traits were measured see the publication methods.

The figures showing the distribution of the phenotypes and the difference between admixture groups are produced with the script [phenotype_analyses_plots.R](code/phenotype_analyses_plots.R).

## Pipeline

### Reads preparation

#### Raw reads download

[01_get_raw_reads.sh](code/01_get_raw_reads.sh) and [rename_reads.sh](code/rename_reads.sh) with reference table [reads_sample_ID.csv](data/reads_sample_ID.csv).

#### Quality control and trimming

[02a_trim_rawreads.sh](code/02a_trim_rawreads.sh) performs quality control with fastqc, then trims the reads with parameters:
` LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 `. It then performs a second quality control.

[02b_fqc_summary.sh](code/02b_fqc_summary.sh) summarises the output of fastqc using a python script: [fastqc_parser.py](code/fastqc_parser.py).
The read counts for raw and trimmed reads is in [raw_reads_count.csv](data/raw_reads_count.csv).


### Alignment

Genome index and dictionary done with script [03_index_dictionary.sh](code/03_index_dictionary.sh).

Alignment done with BWA, script [04_align.sh](code/04_align.sh), with BWA-MEM default parameters.

Duplicated reads marking done with [05a_markdup_metrics_coverage.sh](code/05a_markdup_metrics_coverage.sh), which also calculate some stats about the read mapped and extracts the coordinates of the regions with coverage higher than 100 reads. These regions were then excluded from the variant calling together with the repetitive regions. Repetitive regions identified with repeatMasker, including the modeled regions on the reference genome and [the TREP complete nucleotide database release 19](https://trep-db.uzh.ch/). Step performed in script [05b_mark_repetitive_regions.sh](code/05b_mark_repetitive_regions.sh). The regions to exclude are merged with the annotated repetitive regions with an R script [05c_regions_to_mask.R](code/05c_regions_to_mask.R).

Some numbers of mapped reads, duplicated reads, genome and gene space coverage are listed in [alignment_metrics.csv](data/alignment_metrics.csv).

| Metric | Value | 
| --------------- | --------------- |
| Average coverage genome space | 5.20 |
| Average coverage gene space | 6.20 | 
| Proportion genome covered 1X | 0.79 |
| Proportion genome covered 5X | 0.46 |
| Proportion gene space covered 1X | 0.87|
| Proportion gene space covered 5X | 0.56 |


### Variant calling

Performed with GATK v 4.1.3.0, following best practices. Since Petunia does not have a high quality reference database of variants, we applied filters on the called variants as suggested by GATK. To check that the quality filters suggested in best practices were suitable, we plot the distribution of the quality values of the variants. 

Variant calling performed with script [06_call_variants.sh](code/06_call_variants.sh). We then combined the g.vcf files with script [07_combine_gvcfs.sh](code/07_combine_gvcfs.sh). The genotype calling and variant quality evaluation is done with script [08_genotype_gvcfs_quality.sh](code/08_genotype_gvcfs_quality.sh). In this script we call an R script to plot the distribution of the quality values to verify that the hard filters proposed by GATK are appropriate. This is done with [plot_vcfq_distribution.R](code/plot_vcfq_distribution.R) to make the plots which can be seen in [indel_quality.pdf](data/indel_quality.pdf) and [snp_quality.pdf](data/snp_quality.pdf).

#### Variant filtering

The set of variants is filtered with script [09_filtervcf_quality.sh](code/09_filtervcf_quality.sh).

In this script we filter the set of variants obtained in order to select only those passing quality filters. We then select only variants that are biallelic, and have a call rate higher than 90% across samples, and we exclude the INDELs. We then use ANGSD to calculate the minor allele frequency (MAF) out of the genotype likelihoods, and select only the variants with minor allele frequency > 0.05.

The non-filtered set of variants includes 15 825 083 positions, the quality filtered variants include 13 765 369 positions, the biallelic set includes 12 650 209, the call rate > 90%, only INDELs includes 7 454 530 positions, the MAF > 0.05 includes 4 278 736 positions.

#### Convert vcf to genotype likelihood file (beagle format)

To use the genotype likelihoods in some analyses we had to have them in Beagle format, so I wrote a python script [gatkPLtobeagleGL.py](code/gatkPLtobeagleGL.py) to convert the PL field in the GATK vcf file into genotype likelihoods. Note that missing genoptypes are encoded as 0.33,0.33,0.33 in the output. The conversion is done in script [10_PLtoGL.sh](code/10_PLtoGL.sh).

### Admixture analysis

Is performed with NGSadmix, provided in ANGSD, with script [11a_ngsadmix.sh](code/11a_ngsadmix.sh). I test K from 1 to 8 and each K is run 10 times to then select the best likelihood. The output files `.log` of ngsadmix are parsed to collect the likelihood of each run with the script [ngsadmix_outparser.py](code/ngsadmix_outparser.py):

```
python code/ngsadmix_outparser.py -p data/raw/admixture
```

The results are then analysed in R. with script [11b_ngsadmix.R](code/11b_ngsadmix.R).

### Genomic PCA analysis


### Genetic architecture prediction and GWAS

These are performed with the software [GEMMA](https://github.com/genetics-statistics/GEMMA).

The analyses require phenotype and genotype data as input. To prepare the genotype data I used script [13a_gwas_genotype_formatting.sh]
(code/13a_gwas_genotype_formatting.sh), and to prepare the phenotype I used an R script, [13b_gwas_phenotype_formatting.sh](code/13b_gwas_phenotype_formatting.Rmd). In the latter I also tested the normality of the phenotype distributions.

In the R markdown script [13b](code/13b_gwas_phenotype_formatting.Rmd) I also analysed the phenotype data. I plotted their distributions, and checked their pairwise correlations and correlation to the genomic PCA. The markdown output in html is available in [13b_gwas_phenotype_formatting.html](13b_gwas_phenotype_formatting.html).

I then use the genotype [hardfiltered_biallelic_cr09_mm005.bimbam.gz](data/hardfiltered_biallelic_cr09_mm005.bimbam.gz) and phenotype [pheno_gwas.bimbam](data/pheno_gwas.bimbam) input files to perform the association analyses and the estimation of the genetic architecture. This is done with script [13c_gwas_bslmm.sh](code/13c_gwas_bslmm.sh).

Briefly, we first calculate the kinship matrix with GEMMA, using option `-gk 1` which calculates a centred relatedness matrix.
We then do the association analyses for each of the four phenotypic traits considered (pistil and tube length, anthocyanin and flavonol content), applying first a linear model (LM) with option `-lm 4`, and then a linear mixed model (LMM) `-lmm 4` which includes a correction for relatedness.
All of these analyses are performed with options `-hwe 0.001 -miss 0 -maf 0`.

We then apply the bayesian sparse linear mixed model (BSLMM) to estimate the genetic architecture of each phenotypic trait. This is done with options:

```
-bslmm 1
-hwe 0.001
-maf 0
-miss 0
-w 50000000
-s 200000000
-rpace 1000
-wpace 10000000
```

The output of all analyses is analysed in script [13d_gemma.Rmd](code/13d_gemma.Rmd).


### Divergence scan


## Software versions

- ANGSD/0.933-111-g5859d2b (htslib: 1.11-113-g6038f97)
- bwa/0.7.17
- covtobed/1.2.0 https://github.com/telatin/covtobed
- fastp/0.19.5
- fastqc/0.11.7
- GEMMA/0.98.4
- GenomeAnalysisTK/4.1.3.0
- picard-tools/2.21.8
- R on the computing cluster 3.4.2
- R on the local machine 3.3.3
- RepeatModeler/1.0.11
- samtools/1.10
- sratoolkit/2.10.7
- trimmomatic/0.36
- vcftools/0.1.15

### R libraries

- corrplot
- data.table
- GenomicRanges
- ggplot2
- IRanges
- optparse https://cran.r-project.org/package=optparse
- reshape2
- scales
- 

## References

Bombarely Aureliano, Michel Moser, Avichai Amrad, Manzhu Bao, Laure Bapaume, Cornelius S. Barry, Mattijs Bliek, et al. 2016. ‘Insight into the Evolution of the Solanaceae from the Parental Genomes of Petunia Hybrida’. Nature Plants 2 (6): 16074. https://doi.org/10.1038/nplants.2016.74.

