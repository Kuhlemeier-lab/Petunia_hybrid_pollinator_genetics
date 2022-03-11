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
Alignment done with BWA, script 04_align.sh, parameters:


## Software versions

sratoolkit/2.10.7
fastqc/0.11.7
trimmomatic/0.36
bwa/0.7.17
samtools/1.10
