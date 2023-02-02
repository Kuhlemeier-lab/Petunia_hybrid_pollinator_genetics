# Phenotype descriptive stats, normality, correlation with genetic PCA.
# 23/08/2021
# Last modified 13/01/2023
# Marta Binaghi

wdir <- "hybrids"
setwd(wdir)


# libraries ---------------------------------------------------------------
library(corrplot)
library(reshape2)
library(ggplot2)
library(ggdist) # plot with half violin

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


# import genotype of mybfl -------------------------------------------------------
# only available for contact zone 1
mybfl <- read.table("data/raw/genotype_mybfl_CAPS.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)
pheno_adm_gt <- merge.data.frame(pheno_adm,
                                 mybfl,
                                 by.x = c("parental", "plant"),
                                 by.y = c("parental", "plant"),
                                 all.x = TRUE)

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

# phenotype to phenotype correlation --------------------------------------
# a comprehensive plot with flavonol by anthocyanin, size is exsertion
ggplot(data = pheno,
       aes(x = flav,
           y = antho,
           size = (pist.tube.ratio))) +
  geom_point(alpha = 0.4, col = "gray20") +
  theme_classic() +
  scale_radius(#max_size = 5, trans = "exp"
    #limits = c(0.8, 1.3)
    range = c(0.8, 4))  +
  #theme(legend.position = c(0.87, 0.25)) +
  guides(size = guide_legend(label.position = "bottom")) +
  theme(legend.position="bottom") +
  xlab("Flavonol (AU)") +
  ylab("Anthocyanin (AU)") +
  labs(size = "Exsertion")
ggsave("figures/exploratory/phenotype/scatter_flav_antho_exs.pdf",
       width = 4.4,
       height = 4.7)

# flavonol by antho
ggplot(data = pheno,
       aes(x = flav,
           y = antho)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Flavonol (AU)") +
  ylab("Anthocyanin (AU)")
ggsave("figures/exploratory/phenotype/scatter_flav_antho.pdf",
       width = 2.2,
       height = 4.7/3)
# flavonol by exs
ggplot(data = pheno,
       aes(x = flav,
           y = pist.tube.ratio)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Flavonol (AU)") +
  ylab("Exsertion")
ggsave("figures/exploratory/phenotype/scatter_flav_exs.pdf",
       width = 2.3,
       height = 4.7/3)
# antho by exs
ggplot(data = pheno,
       aes(x = pist.tube.ratio,
           y = antho)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Exsertion") +
  ylab("Anthocyanin (AU)")
ggsave("figures/exploratory/phenotype/scatter_exs_antho.pdf",
       width = 2.2,
       height = 4.7/3)

#correlation across the population
# output file
outfile <- "data/clean/phenotype_correlation.txt"
cat("####  Flavonol - antho   ####",
    file = outfile,
    append = FALSE,
    sep = "\n")
capture.output(cor.test(pheno$flav, pheno$antho),
                 file = outfile,
                 append = TRUE)
cat("####  Flavonol - exsertion   ####",
    file = outfile,
    append = TRUE,
    sep = "\n")
capture.output(cor.test(pheno$flav, pheno$pist.tube.ratio),
               file = outfile,
               append = TRUE)
cat("####  Exsertion - antho   ####",
    file = outfile,
    append = TRUE,
    sep = "\n")
capture.output(cor.test(pheno$pist.tube.ratio, pheno$antho),
               file = outfile,
               append = TRUE)

# phenotype means and SD in admixture groups ----------------------------------------
# phenotypes
phenos <- c("antho", "flav", "pist.tube.ratio")
#adm < 0.10 / >0.90
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
                   median(pheno_adm[pheno_adm$Q2.1 < 0.10, i], na.rm = TRUE))
  mean_adm0 <- c(mean_adm0,
                 mean(pheno_adm[pheno_adm$Q2.1 < 0.10, i], na.rm = TRUE))
  sd_adm0 <- c(sd_adm0,
               sd(pheno_adm[pheno_adm$Q2.1 < 0.10, i], na.rm = TRUE))
  median_adm1 <- c(median_adm1,
                   median(pheno_adm[pheno_adm$Q2.1 > 0.90, i], na.rm = TRUE))
  mean_adm1 <- c(mean_adm1,
                 mean(pheno_adm[pheno_adm$Q2.1 > 0.90, i], na.rm = TRUE))
  sd_adm1 <- c(sd_adm1,
               sd(pheno_adm[pheno_adm$Q2.1 > 0.90, i], na.rm = TRUE))
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
          file = "data/clean/phenotype_means_adm1090.csv", 
          quote = FALSE, 
          row.names = FALSE)


# adm 0.10 0.90
pheno_adm$q0109group <- 0
pheno_adm$q0109group[pheno_adm$Q2.1 > 0.90] <- 1
pheno_adm$q0109group[pheno_adm$Q2.1 < 0.10] <- 2
pheno_adm$q0109group <- as.factor(pheno_adm$q0109group)

# phenotype stats difference between adm groups ---------------------------
#adm 0.10 /0.90
outfile <- "data/clean/phenotype_admGroups1090_statistical_difference.txt"

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
  capture.output(t.test(pheno_adm[pheno_adm$Q2.1 < 0.10, thispheno],
                        pheno_adm[pheno_adm$Q2.1 > 0.90, thispheno],
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
  capture.output(wilcox.test(pheno_adm[pheno_adm$Q2.1 < 0.10, thispheno],
                             pheno_adm[pheno_adm$Q2.1 > 0.90, thispheno],
                             alternative = "two.sided"),
                 file = outfile,
                 append = TRUE)
}



# phenotype in admixture boxplots -----------------------------------------
# adm 0.10/0.90
adm_pure <- rep(NA, times = 70)
adm_pure[pheno_adm$Q2.1 < 0.10] <- "< 0.10"
adm_pure[pheno_adm$Q2.1 > 0.90] <- "> 0.90"
pheno_adm$pure_groups <- adm_pure
pheno_adm$pure_groups <- as.factor(pheno_adm$pure_groups)

# phenotype in admixture violin/boxplot -----------------------------------
pheno_adm$q0109group <- ordered(pheno_adm$q0109group, levels = c("2", "0", "1"))

admixture_col <- c("2" = "#545454",
                   "0" = "#ffc107",
                   "1" = "#e6002e")

admixture_fill <- c("2" = "#545454",
                   "0" = "#ffc107",
                   "1" = "#e6002e")
#antho
ggplot(pheno_adm, 
       aes(x = q0109group,
           y = antho,
           group = q0109group,
           col = q0109group,
           fill = q0109group)) + 
  ggdist::stat_halfeye(adjust = .5,
                       width = .7, 
                       .width = 0, 
                       justification = -.4, 
                       point_colour = NA, alpha = 1) + 
  geom_boxplot(width = .2,
               outlier.shape = NA,
               alpha = 0.2) +
  geom_point(size = 2.8,
             alpha = .8,
             shape = 16,
             position = position_jitter(seed = 1, width = .17)) +
  scale_y_continuous(trans = "log") + 
  theme_classic() +
  scale_color_manual(values = admixture_col) +
  scale_fill_manual(values = admixture_col) +
  theme(legend.position = "none")
ggsave("figures/exploratory/phenotype/boxplot_violin_vert_antho_adm0109.pdf",
       width = 6.5*2/3*1.1,
       height = 6.5)

#flav with genotype of myb-fl
# adm 0.10 0.90 in genotype df
pheno_adm_gt$q0109group <- 0
pheno_adm_gt$q0109group[pheno_adm_gt$Q2.1 > 0.90] <- 1
pheno_adm_gt$q0109group[pheno_adm_gt$Q2.1 < 0.10] <- 2
pheno_adm_gt$q0109group <- as.factor(pheno_adm_gt$q0109group)
pheno_adm_gt$q0109group <- ordered(pheno_adm_gt$q0109group, levels = c("2", "0", "1"))
# code shapes for genotype
pheno_adm_gt$gt_shape <- 0
pheno_adm_gt$gt_shape[pheno_adm_gt$genotype == "e"] <- 15
pheno_adm_gt$gt_shape[pheno_adm_gt$genotype == "h"] <- 16
pheno_adm_gt$gt_shape[pheno_adm_gt$genotype == "a"] <- 17
pheno_adm_gt$gt_shape[is.na(pheno_adm_gt$genotype)] <- 4
pheno_adm_gt$gt_shape <- as.factor(as.character(pheno_adm_gt$gt_shape))

ggplot(pheno_adm_gt, 
       aes(x = q0109group,
           y = flav,
           group = q0109group,
           col = q0109group,
           fill = q0109group)) + 
  ggdist::stat_halfeye(adjust = .5,
                       width = .7, 
                       .width = 0, 
                       justification = -.4, 
                       point_colour = NA, alpha = 1) + 
  geom_boxplot(width = .2,
               outlier.shape = NA,
               alpha = 0.2) +
  geom_point(aes(shape = gt_shape),
             size = 2.8,
             alpha = .8,
             position = position_jitter(seed = 16, width = .17)) +
  theme_classic() +
  scale_color_manual(values = admixture_col) +
  scale_fill_manual(values = admixture_col) +
  scale_shape_manual(values = c(15, 16, 17, 4)) +
  theme(legend.position = "none")
ggsave("figures/exploratory/phenotype/boxplot_violin_vert_flav_adm0109_mybfl.pdf",
       width = 6.5*2/3,
       height = 6.5)
#exs
ggplot(pheno_adm, 
       aes(x = q0109group,
           y = pist.tube.ratio,
           group = q0109group,
           col = q0109group,
           fill = q0109group)) + 
  ggdist::stat_halfeye(adjust = .5,
                       width = .7, 
                       .width = 0, 
                       justification = -.4, 
                       point_colour = NA, alpha = 1) + 
  geom_boxplot(width = .2,
               outlier.shape = NA,
               alpha = 0.2) +
  geom_point(size = 2.8,
             alpha = .8,
             shape = 16,
             position = position_jitter(seed = 1, width = .17)) +
  theme_classic() +
  scale_color_manual(values = admixture_col) +
  scale_fill_manual(values = admixture_col) +
  theme(legend.position = "none") #+
#scale_x_discrete(#limits = factor(c(1.57, 3.8)),
# breaks = seq(1.6, 3.8, by = .2),
# expand = c(-0.701, 2.001)
#)
ggsave("figures/exploratory/phenotype/boxplot_violin_vert_exs_adm0109.pdf",
       width = 6.5*2/3,
       height = 6.5)

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


# phenotype correlation to admixture --------------------------------------
# list the variables that I want to have in the correlation
phenotypes <- c("pist.tube.ratio", "flav", "antho")
outfile <- "data/clean/phenotype_correlation_admixture.txt"
for (thispheno in phenotypes) {
  cat("####    ####",
      file = outfile,
      append = TRUE,
      sep = "\n")
  cat(thispheno,
      file = outfile,
      append = TRUE,
      sep = "\n")
  capture.output(cor.test(pheno_adm$Q2.1,
                          pheno_adm[ , thispheno]),
                 file = outfile,
                 append = TRUE)
}


# plot the phenotype correlation to admixture
# flavonol
ggplot(data = pheno_adm,
       aes(x = Q2.1,
           y = flav)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Admixture proportion") +
  ylab("Flavonol content (AU)")
ggsave("figures/exploratory/phenotype/scatter_admixture_flav.pdf",
       width = 2.2,
       height = 4.7/3)
# anthocyanin
ggplot(data = pheno_adm,
       aes(x = Q2.1,
           y = antho)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Admixture proportion") +
  ylab("Anthocyanin content (AU)")
ggsave("figures/exploratory/phenotype/scatter_admixture_antho.pdf",
       width = 2.2,
       height = 4.7/3)
#exsertion
ggplot(data = pheno_adm,
       aes(x = Q2.1,
           y = pist.tube.ratio)) +
  geom_point(alpha = 0.4, col = "gray20", size = 0.5) +
  theme_classic() +
  xlab("Admixture proportion") +
  ylab("Pistil exsertion")
ggsave("figures/exploratory/phenotype/scatter_admixture_exsertion.pdf",
       width = 2.2,
       height = 4.7/3)
