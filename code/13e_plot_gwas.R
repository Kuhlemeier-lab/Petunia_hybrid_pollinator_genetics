# Plot the GWAS results
# Last modified 14/06/2022

workdir <- "/xxx/hybrids_peaxiINV/"
setwd(workdir)

# libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)

### Tube length ###

# import results  ------------------------------------------------------
gwas <- read.table("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno1_lmm4.assoc.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   as.is = "rs",
                   stringsAsFactors = FALSE)

# extract chromosome names and IDs
gwas$chr.name <- gsub("(Chr\\d)-.*", "\\1", gwas$rs, perl = TRUE)
gwas$chr.name[grepl("Scaffold", gwas$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", gwas$rs, perl = TRUE)[grepl("Scaffold", gwas$chr.name)]

gwas$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_).*", "\\2", gwas$rs, perl = TRUE))
# is it a chromosome or a scaffold
gwas$feature <- gsub("(Chr|Scaffold).*", "\\1", gwas$rs, perl = TRUE)

# extract variant position
gwas$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            gwas$rs, 
                            perl = TRUE))

# make a cumulative position column to concatenate the scaffolds into a single
gwas <- gwas[order(gwas$chr.num, gwas$pos), ]

cum.dist <- 0
gwas$cum.dist <- 0

for (scf in unique(gwas$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    gwas$cum.dist[gwas$chr.name == scf] <- gwas$pos[gwas$chr.name == scf]
  } else {
    gwas$cum.dist[gwas$chr.name == scf] <- cum.dist + gwas$pos[gwas$chr.name == scf]
    cum.dist <- cum.dist + gwas$pos[gwas$chr.name == scf][length(gwas$pos[gwas$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
gwas$chrFacet <- NA
gwas$chrFacet[gwas$feature == "Chr"] <- gwas$chr.name[gwas$feature == "Chr"]
gwas$chrFacet[gwas$feature == "Scaffold"] <- "Scaffolds"

# calculate Bonferroni-corrected P value
# save P values in a vector
gwas$logp <- -log10(gwas$p_lrt)
# calculate Bonferroni
threshold <- 0.05
thresh_value <- -log10( threshold / length(gwas$logp))

# plot ----------------------------------------------------------------
# whole genome
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp)) +
  geom_point(col = "#737373",
    size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/gwas/manhattan_tube.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)
pl <- ggplot(data = gwas,
       aes(x = cum.dist/1000000,
           y = logp,
           col = as.factor(chrFacet))) +
  geom_point(#col = "#737373",
             size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("#ff1e26", 
                                 "#fe941e",
                                 "#eae400", #"#ffff00",
                                 "#06bd00",
                                 "#001a98",
                                 "#760098",
                                 "#f6aab7",
                                 "#55cdfd"))
ggsave("figures/exploratory/gwas/manhattan_tube_rainbow.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)

### Pistil length ###

# import results  ------------------------------------------------------
gwas <- read.table("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno2_lmm4.assoc.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   as.is = "rs",
                   stringsAsFactors = FALSE)

# extract chromosome names and IDs
gwas$chr.name <- gsub("(Chr\\d)-.*", "\\1", gwas$rs, perl = TRUE)
gwas$chr.name[grepl("Scaffold", gwas$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", gwas$rs, perl = TRUE)[grepl("Scaffold", gwas$chr.name)]

gwas$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_).*", "\\2", gwas$rs, perl = TRUE))
# is it a chromosome or a scaffold
gwas$feature <- gsub("(Chr|Scaffold).*", "\\1", gwas$rs, perl = TRUE)

# extract variant position
gwas$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            gwas$rs, 
                            perl = TRUE))

# make a cumulative position column to concatenate the scaffolds into a single
gwas <- gwas[order(gwas$chr.num, gwas$pos), ]

cum.dist <- 0
gwas$cum.dist <- 0

for (scf in unique(gwas$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    gwas$cum.dist[gwas$chr.name == scf] <- gwas$pos[gwas$chr.name == scf]
  } else {
    gwas$cum.dist[gwas$chr.name == scf] <- cum.dist + gwas$pos[gwas$chr.name == scf]
    cum.dist <- cum.dist + gwas$pos[gwas$chr.name == scf][length(gwas$pos[gwas$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
gwas$chrFacet <- NA
gwas$chrFacet[gwas$feature == "Chr"] <- gwas$chr.name[gwas$feature == "Chr"]
gwas$chrFacet[gwas$feature == "Scaffold"] <- "Scaffolds"

# calculate Bonferroni-corrected P value
# save P values in a vector
gwas$logp <- -log10(gwas$p_lrt)
# calculate Bonferroni
threshold <- 0.05
thresh_value <- -log10( threshold / length(gwas$logp))

# plot ----------------------------------------------------------------
# whole genome
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp)) +
  geom_point(col = "#737373",
             size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/gwas/manhattan_pistil.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp,
                 col = as.factor(chrFacet))) +
  geom_point(#col = "#737373",
    size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("#ff1e26", 
                                "#fe941e",
                                "#eae400", #"#ffff00",
                                "#06bd00",
                                "#001a98",
                                "#760098",
                                "#f6aab7",
                                "#55cdfd"))
ggsave("figures/exploratory/gwas/manhattan_pist_rainbow.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)

### Flavonol content ###

# import results  ------------------------------------------------------
gwas <- read.table("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno3_lmm4.assoc.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   as.is = "rs",
                   stringsAsFactors = FALSE)

# extract chromosome names and IDs
gwas$chr.name <- gsub("(Chr\\d)-.*", "\\1", gwas$rs, perl = TRUE)
gwas$chr.name[grepl("Scaffold", gwas$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", gwas$rs, perl = TRUE)[grepl("Scaffold", gwas$chr.name)]

gwas$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_).*", "\\2", gwas$rs, perl = TRUE))
# is it a chromosome or a scaffold
gwas$feature <- gsub("(Chr|Scaffold).*", "\\1", gwas$rs, perl = TRUE)

# extract variant position
gwas$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            gwas$rs, 
                            perl = TRUE))

# make a cumulative position column to concatenate the scaffolds into a single
gwas <- gwas[order(gwas$chr.num, gwas$pos), ]

cum.dist <- 0
gwas$cum.dist <- 0

for (scf in unique(gwas$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    gwas$cum.dist[gwas$chr.name == scf] <- gwas$pos[gwas$chr.name == scf]
  } else {
    gwas$cum.dist[gwas$chr.name == scf] <- cum.dist + gwas$pos[gwas$chr.name == scf]
    cum.dist <- cum.dist + gwas$pos[gwas$chr.name == scf][length(gwas$pos[gwas$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
gwas$chrFacet <- NA
gwas$chrFacet[gwas$feature == "Chr"] <- gwas$chr.name[gwas$feature == "Chr"]
gwas$chrFacet[gwas$feature == "Scaffold"] <- "Scaffolds"

# calculate Bonferroni-corrected P value
# save P values in a vector
gwas$logp <- -log10(gwas$p_lrt)
# calculate Bonferroni
threshold <- 0.05
thresh_value <- -log10( threshold / length(gwas$logp))

# plot ----------------------------------------------------------------
# whole genome
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp)) +
  geom_point(col = "#737373",
             size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/gwas/manhattan_flavonol.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp,
                 col = as.factor(chrFacet))) +
  geom_point(#col = "#737373",
    size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("#ff1e26", 
                                "#fe941e",
                                "#eae400", #"#ffff00",
                                "#06bd00",
                                "#001a98",
                                "#760098",
                                "#f6aab7",
                                "#55cdfd"))
ggsave("figures/exploratory/gwas/manhattan_flavonol_rainbow.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)

### Anthocyanin content ###

# import results  ------------------------------------------------------
gwas <- read.table("data/raw/gwas/output/hardfiltered_biallelic_cr09_mm005_pheno4_lmm4.assoc.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   as.is = "rs",
                   stringsAsFactors = FALSE)

# extract chromosome names and IDs
gwas$chr.name <- gsub("(Chr\\d)-.*", "\\1", gwas$rs, perl = TRUE)
gwas$chr.name[grepl("Scaffold", gwas$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", gwas$rs, perl = TRUE)[grepl("Scaffold", gwas$chr.name)]

gwas$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_).*", "\\2", gwas$rs, perl = TRUE))
# is it a chromosome or a scaffold
gwas$feature <- gsub("(Chr|Scaffold).*", "\\1", gwas$rs, perl = TRUE)

# extract variant position
gwas$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            gwas$rs, 
                            perl = TRUE))

# make a cumulative position column to concatenate the scaffolds into a single
gwas <- gwas[order(gwas$chr.num, gwas$pos), ]

cum.dist <- 0
gwas$cum.dist <- 0

for (scf in unique(gwas$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    gwas$cum.dist[gwas$chr.name == scf] <- gwas$pos[gwas$chr.name == scf]
  } else {
    gwas$cum.dist[gwas$chr.name == scf] <- cum.dist + gwas$pos[gwas$chr.name == scf]
    cum.dist <- cum.dist + gwas$pos[gwas$chr.name == scf][length(gwas$pos[gwas$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
gwas$chrFacet <- NA
gwas$chrFacet[gwas$feature == "Chr"] <- gwas$chr.name[gwas$feature == "Chr"]
gwas$chrFacet[gwas$feature == "Scaffold"] <- "Scaffolds"

# calculate Bonferroni-corrected P value
# save P values in a vector
gwas$logp <- -log10(gwas$p_lrt)
# calculate Bonferroni
threshold <- 0.05
thresh_value <- -log10( threshold / length(gwas$logp))

# plot ----------------------------------------------------------------
# whole genome
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp)) +
  geom_point(col = "#737373",
             size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/gwas/manhattan_antho.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)
pl <- ggplot(data = gwas,
             aes(x = cum.dist/1000000,
                 y = logp,
                 col = as.factor(chrFacet))) +
  geom_point(#col = "#737373",
    size = 0.2) +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  geom_hline(yintercept = thresh_value, lty = 2, col = "black") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("-log10(P)") +
  theme(strip.background = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("#ff1e26", 
                                "#fe941e",
                                "#eae400", #"#ffff00",
                                "#06bd00",
                                "#001a98",
                                "#760098",
                                "#f6aab7",
                                "#55cdfd"))
ggsave("figures/exploratory/gwas/manhattan_antho_rainbow.png",
       plot = pl,
       width = 9,
       height = 2.8,
       dpi = 600)
