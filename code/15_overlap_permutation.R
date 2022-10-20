# overlap of regions in the GWAS with regions with high divergence
# Marta Binaghi
# last modified 16/06/2022


# Gel et al 2016
# https://doi.org/10.1093/bioinformatics/btv562


workdir <- "/xxx/hybrids/"
setwd(workdir)

library(regioneR) 

set.seed(1604)

# define output text file
outfile <- "data/raw/overlap_gwas_divergence_200kbFst.txt"


# import genome info ------------------------------------
# genome info from fasta fai
genome_raw <- read.table("data/genome/Peax403.fasta.fai", sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE)
genome <- getGenome(genome_raw[ , 1:2])

# fst import ---------------------------------------------------------
# import fst results
fst_file <- "data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded_slidingwindow200kb100kb"
fst <- read.table(fst_file,
                  header = TRUE,
                  sep = "\t",
                  dec = ".",
                  stringsAsFactors = FALSE)
row.names(fst) <- 1:length(row.names(fst))
colnames(fst) <- c("chr",
                   "middle.pos",
                   "n.sites",
                   "value")
fst$chr <- as.factor(fst$chr)

# threshold is 95th percentile.
# make a plot with the distribution of the fst values and the threshold line
fst.thresh <- quantile(fst$value, probs = 0.95)
cat(paste0("95th quantile Fst threshold in genome: ",
           fst.thresh),
    file = outfile,
    append = TRUE,
    sep = "\n")
pdf(file = "figures/exploratory/fst/unfolded_histogram_95th_quantile_200kb.pdf",
    width = 6, height = 4.5)
hist((fst$value),
     xlab = "FST")
abline(v = fst.thresh,
       col = "#e6002e")
dev.off()

# mark which positions are above the threshold, add start and stop at side_size
side_size <- 100000
cat(paste0("Fst, size of padding besides significant regions (on each side): ",
           side_size, " bp"),
    file = outfile,
    append = TRUE,
    sep = "\n")
fst$signif <- fst$value >= fst.thresh
selectedFst <- fst[fst$signif, ]
selectedFst$start <- selectedFst$middle.pos - side_size
selectedFst$end <- selectedFst$middle.pos + side_size
selectedFst <- selectedFst[ , c(1, 6, 7)]
fs <- toGRanges(selectedFst)
cat(paste0("N. regions that are significant in Fst: ", length(fs)),
    file = outfile,
    append = TRUE,
    sep = "\n")

# divide fst by chromosome and make a df for each chromosome with the regions
# passing a chromosome-wise 95th threshold

for (i in 1:7) {
  fst.chr <- fst[fst$chr == paste0("Chr", i), ]
  fst.chr.thresh <- quantile(fst.chr$value, probs = 0.95)
  fst.chr$signif <- fst.chr$value >= fst.chr.thresh
  cat(paste0("95th quantile Fst threshold in chromosome ", 
             i, 
             ": ",
             fst.chr.thresh),
      file = outfile,
      append = TRUE,
      sep = "\n")
  selectedFst.chr <- fst.chr[fst.chr$signif, ]
  selectedFst.chr$start <- selectedFst.chr$middle.pos - side_size
  selectedFst.chr$end <- selectedFst.chr$middle.pos + side_size
  selectedFst.chr <- selectedFst.chr[ , c(1, 6, 7)]
  thisfs <- toGRanges(selectedFst.chr)
  assign(paste0("fs.chr", i),
         thisfs)
  cat(paste0("N. regions that are significant in Fst in chr ",
             i,
             ": ", length(thisfs)),
      file = outfile,
      append = TRUE,
      sep = "\n")
  
}
rm(selectedFst.chr, thisfs)


# loop through the different phenotypes

for (pheno in 1:3) {
  # import GWAS lmm --------------------------------------------------------------
  cat(paste0("####  phenotype ", pheno, "  ####"),
      file = outfile,
      append = TRUE,
      sep = "\n")
  lmm <- read.table(paste0("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno",
                           pheno, "_lmm4.assoc.txt"),
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)
  #split chr and pos
  lmm$chr <- gsub("^(Chr|Scaffold_)(\\d+)[_-].*", "\\1\\2", lmm$rs)
  lmm$pos <- as.numeric(gsub("^.*-(\\d+)$", "\\1", lmm$rs))
  
  # Bonferroni corrected threshold and select the sites passing it
  t_lmm <- 0.05
  t_lmm_corr <- t_lmm / nrow(lmm)
  lmm$signif <- lmm$p_lrt <= t_lmm_corr
  
  # make a new df with only the significant regions, with start and end of region
  # padded of side_size bp
  side_size <- 10000
  cat(paste0("GWAS, size of padding besides significant variants (on each side): ",
             side_size, " bp"),
      file = outfile,
      append = TRUE,
      sep = "\n")
  
  selectedLmm <- lmm[lmm$signif, ]
  selectedLmm$start <- selectedLmm$pos - side_size
  selectedLmm$end <- selectedLmm$pos + side_size
  selectedLmm <- selectedLmm[ , c(1, 18, 19)]
  
  # overlap the GWAS results with the Fst over whole genome  ---------------
  gw <- toGRanges(selectedLmm)
  cat(paste0("N. regions that are significant in GWAS: ", length(gw)),
      file = outfile,
      append = TRUE,
      sep = "\n")
  fs <- toGRanges(selectedFst)
  # how many regions overlap in the two sets?
  n.overlap <- numOverlaps(gw, fs)
  # write it to a file
  cat(paste0("Overlaps with genome-wide divergence on N regions: ", n.overlap),
      file = outfile,
      append = TRUE,
      sep = "\n")
  
  if (n.overlap > 0) {
    pt <- overlapPermTest(A = fs,
                          B = gw,
                          ntimes = 1000,
                          genome = genome,
                          non.overlapping = FALSE,
                          verbose = TRUE,
                          alternative = "greater")
    capture.output(summary(pt), 
                   file = outfile,
                   append = TRUE)
    pdf(file = paste0("figures/exploratory/overlap_test_fst200kb_gwas_pheno",
                      pheno, ".pdf"),
        width = 5, height = 5)
    plot(pt)
    dev.off()
    pdf(file = paste0("figures/exploratory/overlap_summary_fst200kb_gwas_pheno",
                      pheno, ".pdf"),
        width = 9, height = 8)
    overlapGraphicalSummary(A = gw,
                            B = fs)
    dev.off()
    
  }
  
  
  # overlap by chromosome ---------------------------------------------------
  # exclude scaffolds because they are too small to do a permutation analysis
  for (chr in 1:7) {
    cat(paste0("Overlaps with chr-wise divergence outliers of chr ", chr),
        file = outfile,
        append = TRUE,
        sep = "\n")
    genome.chr <- genome[seqnames(genome) == paste0("Chr", chr)]
    gw.chr <- toGRanges(selectedLmm)
    gw.chr <- gw.chr[seqnames(gw.chr) == paste0("Chr", chr)]
    if (length(gw.chr) < 1) {
      cat(paste0("No significant regions in GWAS on chr ", chr),
          file = outfile,
          append = TRUE,
          sep = "\n")
      next
    }
    cat(paste0("GWAS significant regions on chr ", chr, ": ", length(gw.chr)),
        file = outfile,
        append = TRUE,
        sep = "\n")
    fs.chr <- get(paste0("fs.chr", chr))
    # how many regions overlap in the two sets?
    n.overlap <- numOverlaps(gw.chr, fs.chr)
    # write it to a file
    cat(paste0("Overlaps with divergence on N regions: ", n.overlap),
        file = outfile,
        append = TRUE,
        sep = "\n")
    if (n.overlap > 0) {
      pt <- overlapPermTest(A = fs.chr,
                            B = gw.chr,
                            ntimes = 1000,
                            genome = genome.chr,
                            non.overlapping = FALSE,
                            verbose = TRUE,
                            alternative = "greater")
      capture.output(summary(pt), 
                     file = outfile,
                     append = TRUE)
      pdf(file = paste0("figures/exploratory/overlap_test_fst200kb_chr", chr, "_gwas_pheno",
                        pheno, ".pdf"),
          width = 5, height = 5)
      plot(pt)
      dev.off()
      pdf(file = paste0("figures/exploratory/overlap_summary_fst200kb_chr", chr,"_gwas_pheno",
                        pheno, ".pdf"),
          width = 9, height = 8)
      overlapGraphicalSummary(A = gw.chr,
                              B = fs.chr)
      dev.off()
      
    }
  }
  
}
