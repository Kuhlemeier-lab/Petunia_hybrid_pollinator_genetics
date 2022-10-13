# Build a coordinate file of regions with coverage exceeding 100 reads
# Natural hybrids project
# Marta Binaghi
# Created 05/05/2021
# Last modified 19/04/2022

setwd("/xxx/hybrids/")

library(ggplot2)
library(IRanges)
library(GenomicRanges)

# The input files were obtained with the tool covtobed https://github.com/telatin/covtobed
# selecting minimum coverage 100.
# These files have coordinated in bed format, meaning that the stop coordinate
# is non-inclusive.

# import regions with high coverage ---------------------------------------

# import each sample into a df
for ( s in 1:70 ) {
  df <- read.table(file = paste0("data/raw/aligned_reads/", s, "_md_mincov100"),
                   sep = "\t",
                   stringsAsFactors = FALSE,
                   header = FALSE)
  colnames(df) <- c("chr", "start", "stop")
  df$length <- df$stop - df$start
  assign(x = paste0("cov", s), value = df)
  rm(df)
}

# merge the intervals so that I have a final table where the regions with high
# coverage in each sample are merged. Ie if sample 1 has region chr1:10-20
# and sample 2 has region chr1:15-25, in the final table I will have chr1:10-25

c1 <- GRanges(
  seqnames = Rle(cov1$chr),
  ranges = IRanges(cov1$start, cov1$stop)) 

for ( i in 2:70 ) {
  covi <- get(paste0("cov", i))
  ci <- GRanges(
    seqnames = Rle(covi$chr),
    ranges = IRanges(covi$start, covi$stop)) 
  c1 <- union(c1, ci)
}


# merge with the regions from repeat masker.
gff <- read.table("data/genome/Peax403.fasta_repeatModeler_TREP_repeats.gff",
                  sep = "\t",
                  header = FALSE,
                  colClasses = c("character", "NULL", "NULL", "integer", "integer",
                                 "NULL", "NULL", "NULL", "NULL"))
colnames(gff) = c("chr", "start", "stop")

gffR <- GRanges(
  seqnames = Rle(gff$chr),
  ranges = IRanges(gff$start, gff$stop)) 
allranges <- union(gffR, c1)

# save the coordinates in a bed file
out <- data.frame(chr = as.character(seqnames(allranges)),
                  start = start(allranges),
                  end = end(allranges))

ind <- gsub("Chr", "", out$chr)
ind <- gsub("Scaffold_(\\d+)__.*", "\\1", ind, perl = TRUE)
ind <- as.integer(ind)

out <- out[order(ind, out$start), ]

write.table(out,
            "data/genome/Peax403_repetitive_and_min100Xcov.bed", 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)



sum(out$end-out$start) / 1500000000
