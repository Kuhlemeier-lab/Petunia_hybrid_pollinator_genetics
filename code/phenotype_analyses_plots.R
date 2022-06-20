# Phenotype descriptive stats, normality plots, correlation with genetic PCA
# 23/08/2021
# Last modified 13/06/2022
# Marta Binaghi

wdir <- "/xxx/hybrids_peaxiINV/"
setwd(wdir)


# import data -------------------------------------------------------------
pheno <- read.table("data/clean/phenotype_sequenced_individuals.csv",
                    sep = ",",
                    header = TRUE,
                    stringsAsFactors = FALSE)
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


# phenotype distribution --------------------------------------------------

pl <- function(){
  qqnorm(pheno_adm$pist, pch = 1, frame = FALSE, main = "Pistil length")
  qqline(pheno_adm$pist, col = "#e6002e", lwd = 2)
  qqnorm((pheno_adm$Dsum), pch = 1, frame = FALSE, main = "Tube length (D1+D2)")
  qqline((pheno_adm$Dsum), col = "#e6002e", lwd = 2)
  qqnorm((pheno_adm$flav), pch = 1, frame = FALSE, main = "Flavonols")
  qqline((pheno_adm$flav), col = "#e6002e", lwd = 2)
  qqnorm((pheno_adm$antho), pch = 1, frame = FALSE, main = "Anthocyanins")
  qqline((pheno_adm$antho), col = "#e6002e", lwd = 2)
}
pdf("figures/exploratory/phenotype/plot_qq_pheno_distribution.pdf",
    width = 6, height = 6)
par(mfrow = c(2,  2))
pl()
dev.off()
pl()


# phenotype normality test ------------------------------------------------
phenos <- c("antho", "flav", "Dsum", "pist")

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
phenos <- c("antho", "flav", "Dsum", "pist")

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
                   median(pheno_adm[pheno_adm$Q2.2 <= 0.25, i], na.rm = TRUE))
  mean_adm0 <- c(mean_adm0,
                 mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, i], na.rm = TRUE))
  sd_adm0 <- c(sd_adm0,
               sd(pheno_adm[pheno_adm$Q2.2 <= 0.25, i], na.rm = TRUE))
  median_adm1 <- c(median_adm1,
                   median(pheno_adm[pheno_adm$Q2.2 >= 0.75, i], na.rm = TRUE))
  mean_adm1 <- c(mean_adm1,
                 mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, i], na.rm = TRUE))
  sd_adm1 <- c(sd_adm1,
               sd(pheno_adm[pheno_adm$Q2.2 >= 0.75, i], na.rm = TRUE))
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
# phenotypes

outfile <- "data/clean/phenotype_admGroups2575_statistical_difference.txt"

# normally distributed:
phenos <- c("Dsum", "pist")
for (thispheno in phenos) {
  cat("####    ####",
      file = outfile,
      append = TRUE,
      sep = "\n")
  cat(thispheno,
      file = outfile,
      append = TRUE,
      sep = "\n")
  capture.output(t.test(pheno_adm[pheno_adm$Q2.2 <= 0.25, thispheno],
                        pheno_adm[pheno_adm$Q2.2 >= 0.75, thispheno],
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
  capture.output(wilcox.test(pheno_adm[pheno_adm$Q2.2 <= 0.25, thispheno],
                             pheno_adm[pheno_adm$Q2.2 >= 0.75, thispheno],
                             alternative = "two.sided"),
                 file = outfile,
                 append = TRUE)
}

# Phenotype distributions --------------------------------------------------
## anthocyanin
pdf(paste0("figures/exploratory/phenotype/density_antho_mean.pdf"),
    width = 3.25,
    height = 3.15)
plot(density(pheno_adm[ , "antho"], na.rm = T), main = "")
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "antho"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "antho"], na.rm = TRUE)),
       col = c( "grey", "#e6002e"))
dev.off()
## flavonols
pdf(paste0("figures/exploratory/phenotype/density_flav_mean.pdf"),
    width = 3.25,
    height = 3.15)
plot(density(pheno_adm[ , "flav"], na.rm = T), main = "")
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "flav"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "flav"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()
## tube length
pdf(paste0("figures/exploratory/phenotype/density_tube_mean.pdf"),
    width = 3.25,
    height = 3.15)
plot(density(pheno_adm[ , "Dsum"], na.rm = T), main = "")
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "Dsum"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "Dsum"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()
## pistil length
pdf(paste0("figures/exploratory/phenotype/density_pist_mean.pdf"),
    width = 3.25,
    height = 3.15)
plot(density(pheno_adm[ , "pist"], na.rm = T), main = "")
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "pist"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "pist"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()

# plot histograms
## anthocyanin
pdf(paste0("figures/exploratory/phenotype/histogram_antho_mean.pdf"),
    width = 3.25,
    height = 3.15)
hist(pheno_adm[ , "antho"],
     breaks = 20)
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "antho"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "antho"], na.rm = TRUE)),
       col = c( "grey", "#e6002e"))
dev.off()
## flavonols
pdf(paste0("figures/exploratory/phenotype/histogram_flav_mean.pdf"),
    width = 3.25,
    height = 3.15)
hist(pheno_adm[ , "flav"],
     breaks = 20)
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "flav"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "flav"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()
## tube length
pdf(paste0("figures/exploratory/phenotype/histogram_tube_mean.pdf"),
    width = 3.25,
    height = 3.15)
hist(pheno_adm[ , "Dsum"],
     breaks = 20)
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "Dsum"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "Dsum"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()
## pistil length
pdf(paste0("figures/exploratory/phenotype/histogram_pist_mean.pdf"),
    width = 3.25,
    height = 3.15)
hist(pheno_adm[ , "pist"],
     breaks = 20)
abline(v = c(mean(pheno_adm[pheno_adm$Q2.2 >= 0.75, "pist"], na.rm = TRUE),
             mean(pheno_adm[pheno_adm$Q2.2 <= 0.25, "pist"], na.rm = TRUE)),
       col = c("grey", "#e6002e"))
dev.off()

# phenotype to admixture scatter plots --------------------------------------------
pdf("figures/exploratory/phenotype/scatter_admixture_flav.pdf",
    width = 4,
    height = 4)
plot(pheno_adm$Q2.1,
     pheno_adm$flav,
     pch = 16,
     cex = 0.7,
     col = "black",
     xlab = "Admixture proportion",
     ylab = "Trait value",
     ylim = c(10, 120))
dev.off()
pdf("figures/exploratory/phenotype/scatter_admixture_antho.pdf",
    width = 4,
    height = 4)
plot(pheno_adm$Q2.1,
     pheno_adm$antho,
     pch = 16,
     cex = 0.7,
     col = "black",
     xlab = "Admixture proportion",
     ylab = "Trait value",
     ylim = c(0, 63))
dev.off()
pdf("figures/exploratory/phenotype/scatter_admixture_pistil.pdf",
    width = 4,
    height = 4)
plot(pheno_adm$Q2.1,
     pheno_adm$pist,
     pch = 16,
     cex = 0.7,
     col = "black",
     xlab = "Admixture proportion",
     ylab = "Trait value",
     ylim = c(3.8, 6.2))
dev.off()
pdf("figures/exploratory/phenotype/scatter_admixture_tube.pdf",
    width = 4,
    height = 4)
plot(pheno_adm$Q2.1,
     pheno_adm$Dsum,
     pch = 16,
     cex = 0.7,
     col = "black",
     xlab = "Admixture proportion",
     ylab = "Trait value",
     ylim = c(3.3, 5.7))
dev.off()

#boxplot
adm_pure <- rep(NA, times = 70)
adm_pure[pheno_adm$Q2.1 <= 0.25] <- "<= 0.25"
adm_pure[pheno_adm$Q2.1 >= 0.75] <- ">= 0.75"
pheno_adm$pure_groups <- adm_pure
pheno_adm$pure_groups <- as.factor(pheno_adm$pure_groups)

pdf("figures/exploratory/phenotype/boxplot_admixture2575_flav.pdf",
    width = 3,
    height = 4)
boxplot(pheno_adm$flav ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value",
        ylim = c(10, 120))
dev.off()
pdf("figures/exploratory/phenotype/boxplot_admixture2575_antho.pdf",
    width = 3,
    height = 4)
boxplot(pheno_adm$antho ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value",
        ylim = c(0, 63))
dev.off()
pdf("figures/exploratory/phenotype/boxplot_admixture2575_pistil.pdf",
    width = 3,
    height = 4)
boxplot(pheno_adm$pist ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value",
        ylim = c(3.8, 6.2))
dev.off()
pdf("figures/exploratory/phenotype/boxplot_admixture2575_tube.pdf",
    width = 3,
    height = 4)
boxplot(pheno_adm$Dsum ~ pheno_adm$pure_groups,
        col = c("#a2a2a2","#e6002e"),
        xlab = "Admixture proportion",
        ylab = "Trait value",
        ylim = c(3.3, 5.7))
dev.off()

