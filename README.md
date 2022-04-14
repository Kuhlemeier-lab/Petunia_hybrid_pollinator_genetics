# Petunia axillaris x P. exserta natural hybrids, genetic basis of a pollination syndrome

## Pipeline

### Reads preparation

#### Raw reads download

01_get_raw_reads.sh and rename_reads.sh with reference table [reads_sample_ID.csv](reads_sample_ID.csv).

#### Quality control and trimming

02_trim_rawreads.sh performs quality control with fastqc, then trims the reads with parameters:
` LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 `

Then performs a second quality control.

02_fqc_summary.sh summarises the output of fastqc using a python script: fastqc_parser.py.
The read counts for raw and trimmed reads is in [raw_reads_count.csv](data/raw_reads_count.csv).



### Alignment

Genome index done with script 03_index.sh.

Alignment done with BWA, script 04_align.sh, with BWA-MEM default parameters.

Duplicated reads marking done with 05a_markdup_metrics_coverage.sh, which also calculate some stats about the read mapped and extracts the coordinates of the regions with coverage higher than 100 reads. These regions will be then excluded from the variant calling. The regions to exclude are merged with the annotated repetitive regions.

Some numbers of mapped reads, duplicated reads, genome and gene space coverage are listed in [alignment_metrics.csv](data/alignment_metrics.csv). On average the aligned reads have 6.20 mean coverage on the gene space and 5.20 on the genome space. The proportion of genome and gene space covered by 1 read is 0.87

| Metric | Value | 
| --------------- | --------------- |
| Average coverage genome space | 5.20 |
| Average coverage gene space | 6.20 | 
| Proportion genome covered 1X | 0.79 |
| Proportion genome covered 5X | 0.46 |
| Proportion gene space covered 1X | 0.87|
| Proportion gene space covered 5X | 0.56 |






## Software versions

- bwa/0.7.17
- fastqc/0.11.7
- samtools/1.10
- sratoolkit/2.10.7
- trimmomatic/0.36


