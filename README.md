# Petunia axillaris x P. exserta natural hybrids, genetic basis of a pollination syndrome

## Pipeline

### Reads preparation

#### Raw reads download

01_get_raw_reads.sh and rename_reads.sh with reference table reads_sample_ID.csv.

#### Quality control and trimming

02_trim_rawreads.sh performs quality control with fastqc, then trims the reads with parameters:
` LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 `

Then performs a second quality control.

02_fqc_summary.sh summarises the output of fastqc using a python script: fstqc_parser.py.

### Alignment

Genome index done with script 03_index.sh.

Alignment done with BWA, script 04_align.sh, with BWA-MEM default parameters.

Duplicated reads marking done with 05a_markdup_metrics_coverage.sh, which also calculate some stats about the read mapped and extracts the coordinates of the regions with coverage higher than 100 reads. These regions will be then excluded from the variant calling. The regions to exclude are merged with the annotated repetitive regions.




## Software versions

- bwa/0.7.17
- fastqc/0.11.7
- samtools/1.10
- sratoolkit/2.10.7
- trimmomatic/0.36


