# Phenotype descriptive stats, normality, correlation with genetic PCA.
# 23/08/2021
# Last modified 13/10/2022
# Marta Binaghi

wdir <- "hybrids/"
setwd(wdir)


# libraries ---------------------------------------------------------------
library(corrplot)
library(reshape2)
library(ggplot2)

# import data -------------------------------------------------------------
pheno <- read.table("data/clean/phenotype_sequenced_individuals.csv",
                    sep = ",",
                    header = TRUE,
                    stringsAsFactors = FALSE)
# add pistil-tube ratio
pheno$pist.tube.ratio <- pheno$pist / pheno$Dsum
# import admixture information
adm <- read.table("data/raw/admixture/q2_admixture.csv",
                  sep = ",",
                  header = TRUE,
                  stringsAsFactors = FALSE)
# combine phenotype and admixture
pheno_adm <- merge.data.frame(pheno,
                              adm,
                              by.x = "new_id_short",
                              by.y = "id")


# phenotype normality test ------------------------------------------------
phenos <- c("antho", "flav", "pist.tube.ratio")

# output file
outfile <- "data/clean/phenotype_normality.txt"
for (i in phenos) {
  cat("####    ####",
      file = outfile,
      append = TRUE,
      sep = "\n")
  cat(i,
      file = outfile,
      append = TRUE,
      sep = "\n")
  capture.output(shapiro.test(pheno_adm[ , i]),
                 file = outfile,
                 append = TRUE)
}

# Indeed, the flavonol and anthocyanin are non normally distributed.
# They are rank-transformed.
pheno_adm$flav_rank <- rank(pheno_adm$flav, na.last = "keep")
pheno_adm$antho_rank <- rank(pheno_adm$antho, na.last = "keep")

# phenotype means and SD in admixture groups ----------------------------------------
# phenotypes
phenos <- c("antho", "flav", "pist.tube.ratio")

median_whole <- numeric()
mean_whole <- numeric()
sd_whole <- numeric()
median_adm0 <- numeric()
mean_adm0 <- numeric()
sd_adm0 <- numeric()
median_adm1 <- numeric()
mean_adm1 <- numeric()
sd_adm1 <- numeric()

for (i in phenos) {
  median_whole <- c(median_whole,
                    median(pheno_adm[ , i], na.rm = TRUE))
  mean_whole <- c(mean_whole,
                  mean(pheno_adm[ , i], na.rm = TRUE))
  sd_whole <- c(sd_whole,
                sd(pheno_adm[ , i], na.rm = TRUE))
  median_adm0 <- c(median_adm0,
                   median(pheno_adm[pheno_adm$Q2.1 < 0.25, i], na.rm = TRUE))
  mean_adm0 <- c(mean_adm0,
                 mean(pheno_adm[pheno_adm$Q2.1 < 0.25, i], na.rm = TRUE))
  sd_adm0 <- c(sd_adm0,
               sd(pheno_adm[pheno_adm$Q2.1 < 0.25, i], na.rm = TRUE))
  median_adm1 <- c(median_adm1,
                   median(pheno_adm[pheno_adm$Q2.1 > 0.75, i], na.rm = TRUE))
  mean_adm1 <- c(mean_adm1,
                 mean(pheno_adm[pheno_adm$Q2.1 > 0.75, i], na.rm = TRUE))
  sd_adm1 <- c(sd_adm1,
               sd(pheno_adm[pheno_adm$Q2.1 > 0.75, i], na.rm = TRUE))
}

pheno_means <- data.frame(phenotype = phenos,
                          median_whole = median_whole,
                          mean_whole = mean_whole,
                          sd_whole = sd_whole,
                          median_adm0 = median_adm0,
                          mean_adm0 = mean_adm0,
                          sd_adm0 = sd_adm0,
                          median_adm1 = median_adm1,
                          mean_adm1 = mean_adm1,
                          sd_adm1 = sd_adm1)
write.csv(pheno_means, 
          file = "data/clean/phenotype_means_adm2575.csv", 
          quote = FALSE, 
          row.names = FALSE)

# phenotype stats difference between adm groups ---------------------------

outfile <- "data/clean/phenotype_admGroups2575_statistical_difference.txt"

# normally distributed:
phenos <- c("pist.tube.ratio")
for (thispheno in phenos) {
  cat("####    ####",
      file = outfile,
      append = TRUE,
      sep = "\n")
  cat(thispheno,
      file = outfile,
      append = TRUE,
      sep = "\n")
  capture.output(t.test(pheno_adm[pheno_adm$Q2.1 < 0.25, thispheno],
                        pheno_adm[pheno_adm$Q2.1 > 0.75, thispheno],
                        alternative = "two.sided"),
                 file = outfile,
                 append = TRUE)
}

# non-normally distributed:
phenos <- c("antho", "flav")

for (thispheno in phenos) {
  cat("####    ####",
      file = outfile,
      append = TRUE,
      sep = "\n")
  cat(thispheno,
      file = outfile,
      append = TRUE,
      sep = "\n")
  capture.output(wilcox.test(pheno_adm[pheno_adm$Q2.1 < 0.25, thispheno],
                             pheno_adm[pheno_adm$Q2.1 > 0.75, thispheno],
                             alternative = "two.sided"),
                 file = outfile,
                 append = TRUE)
}

# phenotype in admixture boxplots -----------------------------------------
adm_pure <- rep(NA, times = 70)
adm_pure[pheno_adm$Q2.1 < 0.25] <- "< 0.25"
adm_pure[pheno_adm$Q2.1 > 0.75] <- "> 0.75"
pheno_adm$pure_groups <- adm_pure
pheno_adm$pure_groups <- as.factor(pheno_adm$pure_groups)

pdf("figures/exploratory/phenotype/boxplot_admixture2575_flav_thin.pdf",
    width = 2,
    height = 4)
boxplot(pheno_adm$flav ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value")
dev.off()
pdf("figures/exploratory/phenotype/boxplot_admixture2575_antho_thin.pdf",
    width = 2,
    height = 4)
boxplot(pheno_adm$antho ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value")
dev.off()
pdf("figures/exploratory/phenotype/boxplot_admixture2575_pist.tube.ratio_thin.pdf",
    width = 2,
    height = 4)
boxplot(pheno_adm$pist.tube.ratio ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value")
dev.off()

# Phenotype correlation to principal components ---------------------------
# I test phenotype correlation to the first 10 genomic principal components.

# import PCA data
covmat <- read.table("data/raw/pca/cr09_mm005.cov")
# the eigen vector of the covariance matrix:
e <- eigen(covmat)
# the variance explained by each PC:
PCperc <- round(e$values/sum(e$values) * 100, digits = 2)
PCperc[1:10]

# take first 10 PCs
pcs1.10 <- e$vectors[, 1:10]
# sort pheno df as the PCA
pheno_adm <- pheno_adm[order(pheno_adm$id_gl), ]

# list the variables that I want to have in the correlation matrix
phenotypes <- c("pist.tube.ratio",
                "flav", "antho")

# make one matrix to keep the correlation estimate value
estimates_mt <- matrix(data = 999, nrow = length(phenotypes), ncol = 10)
# and one to keep the corresponding P values
pvals_mt <- matrix(data = 999, nrow = length(phenotypes), ncol = 10)

# loop through the phenotypes and calculate their correlation with each PC
for ( thispheno in phenotypes ) {
  test <- apply(X = pcs1.10, 
                MARGIN = 2, 
                function(x) cor.test(as.numeric(pheno_adm[ , thispheno]), x))
  estimates <- sapply(test, "[[", "estimate")
  pvals <- sapply(test, "[[", "p.value")
  estimates_mt[match(thispheno, phenotypes), ] <- estimates
  pvals_mt[match(thispheno, phenotypes), ] <- pvals
}
# name rows and cols of the matrices
colnames(estimates_mt) <- paste0("PC", 1:10)
rownames(estimates_mt) <- phenotypes
colnames(pvals_mt) <- paste0("PC", 1:10)
rownames(pvals_mt) <- phenotypes

# apply Bonferroni correction to the P values by multiplying the P values
# by the number of comparisons performed
pvals_bonf_mt <- pvals_mt * ( length(phenotypes) * 10 )

# make the matrices into a long df
estimates_mt.long <- melt(estimates_mt)
pvals_mt.long <- melt(pvals_mt)
pvals_bonf_mt.long <- melt(pvals_bonf_mt)
corr_PCpheno <- merge.data.frame(estimates_mt.long, pvals_mt.long, by = c("Var1", "Var2"))
corr_PCpheno <- merge.data.frame(corr_PCpheno, pvals_bonf_mt.long, by = c("Var1", "Var2"))
colnames(corr_PCpheno) <- c("Variable1", "Variable2", "Pearson_r2", "Pvalue", "Bonferroni_Pvalue")

# add a column to code significance of Bonferroni corrected P values
corr_PCpheno$bonf_signif <- "na"
corr_PCpheno$bonf_signif[corr_PCpheno$Bonferroni_Pvalue < 0.01] <- "**"
corr_PCpheno$bonf_signif[corr_PCpheno$Bonferroni_Pvalue < 0.05 &
                           corr_PCpheno$Bonferroni_Pvalue >= 0.01 ] <- "*"
corr_PCpheno$bonf_signif[corr_PCpheno$Bonferroni_Pvalue >= 0.05] <- "n.s."

pl <- ggplot(data = corr_PCpheno,
             aes(Variable1, Variable2)) +
  geom_tile(aes(fill = Pearson_r2),
            colour = "white") +        # adds the thin white line bw cells
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  geom_text(aes(label = bonf_signif)) +
  xlab("Phenotypic trait") +
  ylab("PC component") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave("figures/exploratory/phenotype/correlogram_PCTopheno.png",
       plot = pl,
       width = 8, height = 5)
pl

# save the results to a readable csv. These are the results reported in
# Additional file 1, table S2
write.csv(x = corr_PCpheno, 
          file = "data/raw/phenotype/correlation_PCPheno.csv", 
          quote = FALSE, 
          row.names = FALSE)

rm(estimates_mt, 
   estimates_mt.long, 
   pvals_bonf_mt, 
   pvals_bonf_mt.long, 
   pvals_mt, 
   pvals_mt.long, 
   test,
   thispheno,
   pvals,
   estimates)
