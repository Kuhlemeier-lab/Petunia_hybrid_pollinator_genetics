# overlap of regions in the GWAS with regions with high divergence
# Marta Binaghi
# last modified 21/12/2022


# Gel et al 2016
# https://doi.org/10.1093/bioinformatics/btv562


workdir <- "hybrids"
setwd(workdir)

library(regioneR) # note that on R v3 the function seems to have a bug so I 
# corrected it. See at the bottom of this script in case.

set.seed(1604)

# define output text file
outfile <- "data/raw/overlap_gwas_divergence_200kbFst.txt"


# import genome info ------------------------------------
# genome info from fasta fai
genome_raw <- read.table("data/genome/Peax403.fasta.fai", sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE)
genome <- getGenome(genome_raw[ , 1:2])

# adm 0.10 0.90 -----------------------------------------------------------
# fst import ---------------------------------------------------------
outfile <- "data/raw/overlap_gwas_divergence0109_200kbFst.txt"
# import fst results
fst_file <- "data/raw/fst/hf_ba_cr09_mm005_admixture0901_unfolded_slidingwindow200kb100kb"
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
pdf(file = "figures/exploratory/fst/unfolded0109_histogram_95th_quantile_200kb.pdf",
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
for (pheno in 1:5) {
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
  # save the regions to a tab
  write.csv(selectedLmm,
            file = paste0("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno",
                          pheno, "_lmm4.assoc_sitesAbove",
                          t_lmm,
                          ".csv"),
            row.names = FALSE,
            quote = FALSE)
  
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
    pdf(file = paste0("figures/exploratory/overlap_test_fst200kb_adm0109_gwas_pheno",
                      pheno, ".pdf"),
        width = 5, height = 5)
    plot(pt)
    dev.off()
    pdf(file = paste0("figures/exploratory/overlap_summary_fst200kb_adm0109_gwas_pheno",
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


# test overlap with selection of 10 highest sites in the GWAS -------------
top_n <- 10
outfile <- paste0("data/raw/overlap_gwasTop", top_n, "sites_divergence0109_200kbFst.txt")

for (pheno in 3:5) {
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
  
  # define P-value threshold for significance based on the top_n number
  # i.e. if top_n is 10, we want the lowest 10 P-values, so we calculate the threshold 
  # based on the P-value of the 11th lowest P.
  thresh_p <- sort(lmm$p_lrt)[top_n+1]
  cat(paste0("Considering the ", top_n, " sites with lowest P-value.
These correspond to a P value threshold of ", thresh_p),
file = outfile,
append = TRUE,
sep = "\n")
  lmm$signif <- lmm$p_lrt < thresh_p
  
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
  # add the fst of the window that overlaps the GWAS site
  selectedLmm$middle <- ((selectedLmm$end - selectedLmm$start) / 2 ) + 
    selectedLmm$start
  
  selectedLmm$fst.mid.pos <- 0
  selectedLmm$fst.value <- 0
  
  for (i in 1:nrow(selectedLmm)) {
    chri <- selectedLmm$chr[i]
    midi <- selectedLmm$middle[i]
    fsti <- fst.chr$value[which(abs(fst.chr$middle.pos - midi) == min(abs(fst.chr$middle.pos - midi)))]
    fsti.pos <- fst.chr$middle.pos[which(abs(fst.chr$middle.pos - midi) == min(abs(fst.chr$middle.pos - midi)))]
    selectedLmm$fst.mid.pos[i] <- fsti.pos
    selectedLmm$fst.value[i] <- fsti
  }
  
  
  # save the regions to a tab
  write.csv(selectedLmm,
            file = paste0("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno",
                          pheno, "_lmm4.assoc_top", top_n, "sites.csv"),
            row.names = FALSE,
            quote = FALSE)
  
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
    pdf(file = paste0("figures/exploratory/overlap_test_fst200kb_adm0109_gwas_pheno",
                      pheno, ".pdf"),
        width = 5, height = 5)
    plot(pt)
    dev.off()
    pdf(file = paste0("figures/exploratory/overlap_summary_fst200kb_adm0109_gwas_pheno",
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


# session info ------------------------------------------------------------

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Italian_Italy.1252  LC_CTYPE=Italian_Italy.1252    LC_MONETARY=Italian_Italy.1252
# [4] LC_NUMERIC=C                   LC_TIME=Italian_Italy.1252    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] regioneR_1.26.1      GenomicRanges_1.46.1 GenomeInfoDb_1.30.1  IRanges_2.28.0       S4Vectors_0.32.3    
# [6] BiocGenerics_0.40.0 
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.1.2              restfulr_0.0.13             XVector_0.34.0              MatrixGenerics_1.6.0       
# [5] bitops_1.0-7                tools_4.1.2                 zlibbioc_1.40.0             memoise_2.0.1              
# [9] BSgenome_1.62.0             lattice_0.20-45             rlang_1.0.1                 Matrix_1.4-0               
# [13] DelayedArray_0.20.0         cli_3.1.1                   rstudioapi_0.13             yaml_2.2.2                 
# [17] parallel_4.1.2              fastmap_1.1.0               GenomeInfoDbData_1.2.7      rtracklayer_1.54.0         
# [21] Biostrings_2.62.0           grid_4.1.2                  Biobase_2.54.0              XML_3.99-0.8               
# [25] BiocParallel_1.28.3         Rsamtools_2.10.0            matrixStats_0.61.0          GenomicAlignments_1.30.0   
# [29] SummarizedExperiment_1.24.0 RCurl_1.98-1.6              cachem_1.0.6                crayon_1.5.0               
# [33] rjson_0.2.21                BiocIO_1.4.0 

# permTest function Rv3 ---------------------------------------------------

# functions ----------------------------------------------------------------
# # redefine permTest to solve bug with dimension of random.evaluate[, i]
# mypermTest <- function (A, ntimes = 100, randomize.function, evaluate.function, 
#                     alternative = "auto", min.parallel = 1000, force.parallel = NULL, 
#                     randomize.function.name = NULL, evaluate.function.name = NULL, 
#                     verbose = FALSE, ...) 
# {
#   alternative <- match.arg(alternative, c("less", "greater", 
#                                           "auto"))
#   if (!hasArg(A)) 
#     stop("A is missing")
#   if (!is.numeric(ntimes)) 
#     stop("ntime must be numeric")
#   if (!hasArg(randomize.function)) 
#     stop("randomize.function is missing")
#   if (!is.function(randomize.function)) 
#     stop("randomize.function must be a function")
#   if (!hasArg(evaluate.function)) 
#     stop("evaluate.function is missing")
#   if (!(is.function(evaluate.function) | is.list(evaluate.function))) 
#     stop("evaluate.function must be a function")
#   if (!is.numeric(min.parallel)) 
#     stop("min.parallel must be numeric")
#   if (ntimes < 100) 
#     print(paste0("Note: The minimum p-value with only ", 
#                  ntimes, " permutations is ", 1/(ntimes + 1), ". You should consider increasing the number of permutations."))
#   A <- toGRanges(A)
#   if (!is.null(force.parallel)) {
#     doParallel <- force.parallel
#   }
#   else {
#     doParallel <- (length(A) * ntimes > min.parallel)
#   }
#   if (!is.list(evaluate.function)) {
#     if (is.null(evaluate.function.name)) {
#       evaluate.function.name <- as.character(match.call()["evaluate.function"])
#     }
#     ef <- list()
#     ef[[evaluate.function.name]] <- evaluate.function
#     evaluate.function <- ef
#   }
#   else {
#     if (!is.null(evaluate.function.name)) {
#       names(evaluate.function) <- evaluate.function.name
#     }
#     else {
#       if (is.null(names(evaluate.function))) {
#         names(evaluate.function) <- paste0("Function", 
#                                            c(1:length(evaluate.function)))
#       }
#     }
#   }
#   if (is.null(randomize.function.name)) {
#     randomize.function.name <- match.call()["randomize.function"]
#   }
#   original.evaluate <- sapply(c(1:length(evaluate.function)), 
#                               function(i, ...) {
#                                 return(evaluate.function[[i]](A, ...))
#                               }, ...)
#   if (!is.numeric(original.evaluate)) {
#     stop(paste0("The evaluation function must return a numeric value but it returned an object of class ", 
#                 class(original.evaluate)))
#   }
#   if (verbose) {
#     pb <- txtProgressBar(min = 0, max = ntimes, style = 3)
#     setTxtProgressBar(pb, 0)
#   }
#   randomize_and_evaluate <- function(foo, ...) {
#     randomA <- randomize.function(A, ...)
#     if (verbose) {
#       setTxtProgressBar(pb, foo)
#     }
#     rand.evaluate <- sapply(c(1:length(evaluate.function)), 
#                             function(i, ...) {
#                               return(evaluate.function[[i]](randomA, ...))
#                             }, ...)
#     return(rand.evaluate)
#   }
#   if (doParallel) {
#     if (verbose) {
#       random.evaluate <- numeric()
#       chunk.size <- max(round(ntimes/100 + 1), 10)
#       e <- 0
#       done <- FALSE
#       while (!done) {
#         s <- e + 1
#         e <- s + chunk.size
#         if (e >= ntimes) {
#           e <- ntimes
#           done <- TRUE
#         }
#         random.evaluate <- c(random.evaluate, do.call(rbind, 
#                                                       mclapply(c(s:e), randomize_and_evaluate, ...)))
#         setTxtProgressBar(pb, e)
#       }
#     }
#     else {
#       random.evaluate <- do.call(rbind, mclapply(c(1:ntimes), 
#                                                  randomize_and_evaluate, ...))
#     }
#   }
#   else {
#     random.evaluate <- do.call(rbind, lapply(c(1:ntimes), 
#                                              randomize_and_evaluate, ...))
#   }
#   results <- list()
#   for (i in c(1:length(evaluate.function))) {
#     func.name <- names(evaluate.function)[i]
#     orig.ev <- original.evaluate[i]
#     # marta's edit
#     if ( is.vector(random.evaluate) ) {
#       rand.ev <- random.evaluate
#     } else {
#       rand.ev <- random.evaluate[, i]
#     }
#     num.nas <- length(which(is.na(rand.ev)))
#     num.valid.values <- ntimes - num.nas
#     if (num.valid.values < ntimes) {
#       if (num.valid.values > 0) {
#         warning(paste0(num.nas, " iterations returned NA or NaN. Only ", 
#                        , " iterations have been used to compute the p-value."))
#       }
#       else {
#         warning(paste0("All ", num.nas, " iterations returned NA or NaN. No valid values returned. It is not possible to compute the p-value nor z-score."))
#       }
#     }
#     if (num.valid.values > 0) {
#       if (alternative == "auto") {
#         alt <- ifelse(orig.ev < mean(rand.ev, na.rm = TRUE), 
#                       "less", "greater")
#       }
#       else {
#         alt <- alternative
#       }
#       if (alt == "less") {
#         pval <- (sum(orig.ev >= rand.ev, na.rm = TRUE) + 
#                    1)/(num.valid.values + 1)
#       }
#       else {
#         pval <- (sum(orig.ev <= rand.ev, na.rm = TRUE) + 
#                    1)/(num.valid.values + 1)
#       }
#       if (alternative == "greater" & orig.ev < mean(rand.ev, 
#                                                     na.rm = TRUE)) 
#         message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
#       if (alternative == "less" & orig.ev > mean(rand.ev, 
#                                                  na.rm = TRUE)) 
#         message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
#       if (orig.ev == 0 & all(rand.ev == 0)) {
#         warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
#         pval <- 1
#         zscore <- NA
#       }
#       else {
#         zscore <- round((orig.ev - mean(rand.ev, na.rm = TRUE))/sd(rand.ev, 
#                                                                    na.rm = TRUE), 4)
#       }
#     }
#     else {
#       pval <- NA
#       zscore <- NA
#       alt <- alternative
#     }
#     if (!is.na(pval)) {
#       if (!is.finite(zscore)) {
#         warning(paste0("All permuted values are equal to ", 
#                        rand.ev[1], ". Z-score is infinite."))
#       }
#     }
#     res <- list(pval = pval, ntimes = ntimes, alternative = alt, 
#                 observed = orig.ev, permuted = rand.ev, zscore = zscore, 
#                 evaluate.function = evaluate.function[[i]], evaluate.function.name = func.name, 
#                 randomize.function = randomize.function, randomize.function.name = randomize.function.name)
#     class(res) <- "permTestResults"
#     results[[func.name]] <- res
#   }
#   class(results) <- "permTestResultsList"
#   return(results)
# }
# environment(mypermTest) <- asNamespace('regioneR')
# assignInNamespace("permTest", mypermTest, ns = "regioneR")
