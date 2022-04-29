# Petunia axillaris x P. exserta natural hybrids, genetic basis of a pollination syndrome

Bioinformatic analyses performed in the experiments for the publication about the genetic basis of pollination syndromes in Petunia axillaris and P. exserta.

Project conducted in the [Plant genetics and development group](https://www.ips.unibe.ch/research/deve/index_eng.html) at the Institute of Plant Sciences of the University of Bern.

Authors: Marta Binaghi, Marius Roesti, Therese Mandel, Korinna Esfeld, Loreta B. Freitas, Cris Kuhlemeier

Author of this page: Marta Binaghi

## DNA sequencing

Library preparation and sequencing were performed by the [Next Generation Sequencing platform of the University of Bern](https://www.ngs.unibe.ch/) in two batches, one in 2016 and one in 2018. In both batches, the DNA was amplified with illustraTM GenomiPhiTM V2 DNA Amplification Kit. Library preparation followed the TruSeq DNA PCR-free protocol. Sequencing was performed to obtain 150 bp long, paired-end reads, for an estimated coverage of 4-5 ×, on a genome size of 1.2 Gb (Bombarely et al., 2016). The 2016 batch was sequenced on two lanes of an Illumina HiSeq 3000. The 2018 batch was sequenced on two lanes (one chip) of the Illumina NovaSeq.

Raw reads were uploaded to NCBI SRA under [BioProjects PRJNA522653 (2016 batch)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA522653) and [PRJNA706535 (2018 batch)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA706535).


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

Variant calling performed with script [06_call_variants.sh](code/06_call_variants.sh). We then combined the g.vcf files with script [07_combine_gvcfs.sh](code/07_combine_gvcfs.sh). The genotype calling and variant quality evaluation is done with script [08_genotype_gvcfs_quality.sh](08_genotype_gvcfs_quality.sh). In this script we call an R script to plot the distribution of the quality values to verify that the hard filters proposed by GATK are appropriate. This is done with [plot_vcfq_distribution.R](code/plot_vcfq_distribution.R) to make the plots which can be seen in [indel_quality.pdf](data/indel_quality.pdf) and [snp_quality.pdf](data/snp_quality.pdf).

#### Variant filtering

The set of variants is filtered with script [09_filtervcf_quality.sh](code/09_filtervcf_quality.sh).

In this script we filter the set of variants obtained in order to select only those passing quality filters. We then select only variants that are biallelic, and have a call rate higher than 90% across samples, and we exclude the INDELs. We then use ANGSD to calculate the minor allele frequency (MAF) out of the genotype likelihoods, and select only the variants with minor allele frequency > 0.05.

The non-filtered set of variants includes 15 825 083 positions, the quality filtered variants include 13 765 369 positions, the biallelic set includes 12 650 209, the call rate > 90%, only INDELs includes 7 454 530 positions, the MAF > 0.05 includes 4 278 736 positions.

#### Convert vcf to genotype likelihood file (beagle format)

To use the genotype likelihoods in some analyses we had to have them in Beagle format, so I wrote a python script [gatkPLtobeagleGL.py](code/gatkPLtobeagleGL.py) to convert the PL field in the GATK vcf file into genotype likelihoods. Note that missing genoptypes are encoded as 0.33,0.33,0.33 in the output. The conversion is done in script [10_PLtoGL.sh](code/10_PLtoGL.sh).

## Software versions

- bwa/0.7.17
- covtobed/1.2.0 https://github.com/telatin/covtobed
- fastp/0.19.5
- fastqc/0.11.7
- GenomeAnalysisTK/4.1.3.0
- picard-tools/2.21.8
- R on the computing cluster 3.4.2
- R on the local machine 3.3.3
- RepeatModeler/1.0.11
- samtools/1.10
- sratoolkit/2.10.7
- trimmomatic/0.36

### R libraries

- optparse https://cran.r-project.org/package=optparse
- 

## References

Bombarely, Aureliano, Michel Moser, Avichai Amrad, Manzhu Bao, Laure Bapaume, Cornelius S. Barry, Mattijs Bliek, et al. 2016. ‘Insight into the Evolution of the Solanaceae from the Parental Genomes of Petunia Hybrida’. Nature Plants 2 (6): 16074. https://doi.org/10.1038/nplants.2016.74.

