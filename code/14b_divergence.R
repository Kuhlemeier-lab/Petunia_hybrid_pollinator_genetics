# this script uses as input the files produced with realSFS
# Last modified 16/06/2022

workdir <- "/xxx/hybrids_peaxiINV/"
setwd(workdir)

# libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)

# import fst unfolded ------------------------------------------------------
fst <- read.table("data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded_slidingwindow50kb25kb",
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

# extract chromosome names and IDs
fst$chr.name <- gsub("(Chr\\d).*", "\\1", fst$chr, perl = TRUE)
fst$chr.name[grepl("Scaffold", fst$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", fst$chr, perl = TRUE)[grepl("Scaffold", fst$chr.name)]

fst$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_)?.*", "\\2", fst$chr.name, perl = TRUE))
# is it a chromosome or a scaffold
fst$feature <- gsub("(Chr|Scaffold).*", "\\1", fst$chr.name, perl = TRUE)

# make a cumulative position column to concatenate the scaffolds into a single
fst <- fst[order(fst$chr.num, fst$middle.pos), ]

cum.dist <- 0
fst$cum.dist <- 0

for (scf in unique(fst$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    fst$cum.dist[fst$chr.name == scf] <- fst$middle.pos[fst$chr.name == scf]
  } else {
    fst$cum.dist[fst$chr.name == scf] <- cum.dist + fst$middle.pos[fst$chr.name == scf]
    cum.dist <- cum.dist + fst$middle.pos[fst$chr.name == scf][length(fst$middle.pos[fst$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
fst$chrFacet <- NA
fst$chrFacet[fst$feature == "Chr"] <- fst$chr.name[fst$feature == "Chr"]
fst$chrFacet[fst$feature == "Scaffold"] <- "Scaffolds"

fst$mb.pos <- fst$cum.dist/1000000

# plot ----------------------------------------------------------------
# whole genome
ggplot(data = fst,
       aes(x = cum.dist/1000000,
           y = value)) +
  geom_point(col = "#737373",
             size = 0.2) +
  geom_smooth(col = "#e6002e") +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("Fst value") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/fst/fst_admixturePops2575_unfolded_wind50kstep25k.png",
       width = 8,
       height = 3,
       dpi = 600)

# import fst folded ------------------------------------------------------
fst <- read.table("data/raw/fst/hf_ba_cr09_mm005_admixture_folded_slidingwindow50kb25kb",
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

# extract chromosome names and IDs
fst$chr.name <- gsub("(Chr\\d).*", "\\1", fst$chr, perl = TRUE)
fst$chr.name[grepl("Scaffold", fst$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", fst$chr, perl = TRUE)[grepl("Scaffold", fst$chr.name)]

fst$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_)?.*", "\\2", fst$chr.name, perl = TRUE))
# is it a chromosome or a scaffold
fst$feature <- gsub("(Chr|Scaffold).*", "\\1", fst$chr.name, perl = TRUE)

# make a cumulative position column to concatenate the scaffolds into a single
fst <- fst[order(fst$chr.num, fst$middle.pos), ]

cum.dist <- 0
fst$cum.dist <- 0

for (scf in unique(fst$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    fst$cum.dist[fst$chr.name == scf] <- fst$middle.pos[fst$chr.name == scf]
  } else {
    fst$cum.dist[fst$chr.name == scf] <- cum.dist + fst$middle.pos[fst$chr.name == scf]
    cum.dist <- cum.dist + fst$middle.pos[fst$chr.name == scf][length(fst$middle.pos[fst$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
fst$chrFacet <- NA
fst$chrFacet[fst$feature == "Chr"] <- fst$chr.name[fst$feature == "Chr"]
fst$chrFacet[fst$feature == "Scaffold"] <- "Scaffolds"

fst$mb.pos <- fst$cum.dist/1000000

# plot ----------------------------------------------------------------
# whole genome
ggplot(data = fst,
       aes(x = cum.dist/1000000,
           y = value)) +
  geom_point(col = "#737373",
             size = 0.2) +
  geom_smooth(col = "#e6002e") +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("Fst value") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/fst/fst_admixturePops2575_folded_wind50kstep25k.png",
       width = 8,
       height = 3,
       dpi = 600)


# 100 kb window -----------------------------------------------------------
fst <- read.table("data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded_slidingwindow100kb50kb",
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

# extract chromosome names and IDs
fst$chr.name <- gsub("(Chr\\d).*", "\\1", fst$chr, perl = TRUE)
fst$chr.name[grepl("Scaffold", fst$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", fst$chr, perl = TRUE)[grepl("Scaffold", fst$chr.name)]

fst$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_)?.*", "\\2", fst$chr.name, perl = TRUE))
# is it a chromosome or a scaffold
fst$feature <- gsub("(Chr|Scaffold).*", "\\1", fst$chr.name, perl = TRUE)

# make a cumulative position column to concatenate the scaffolds into a single
fst <- fst[order(fst$chr.num, fst$middle.pos), ]

cum.dist <- 0
fst$cum.dist <- 0

for (scf in unique(fst$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    fst$cum.dist[fst$chr.name == scf] <- fst$middle.pos[fst$chr.name == scf]
  } else {
    fst$cum.dist[fst$chr.name == scf] <- cum.dist + fst$middle.pos[fst$chr.name == scf]
    cum.dist <- cum.dist + fst$middle.pos[fst$chr.name == scf][length(fst$middle.pos[fst$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
fst$chrFacet <- NA
fst$chrFacet[fst$feature == "Chr"] <- fst$chr.name[fst$feature == "Chr"]
fst$chrFacet[fst$feature == "Scaffold"] <- "Scaffolds"

fst$mb.pos <- fst$cum.dist/1000000

# plot ----------------------------------------------------------------
# whole genome
ggplot(data = fst,
       aes(x = cum.dist/1000000,
           y = value)) +
  geom_point(col = "#737373",
             size = 0.2) +
  geom_smooth(col = "#e6002e") +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("Fst value") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/fst/fst_admixturePops2575_unfolded_wind100kstep50k.png",
       width = 8,
       height = 3,
       dpi = 600)

# 200 kb window -----------------------------------------------------------
fst <- read.table("data/raw/fst/hf_ba_cr09_mm005_admixture_unfolded_slidingwindow200kb100kb",
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

# extract chromosome names and IDs
fst$chr.name <- gsub("(Chr\\d).*", "\\1", fst$chr, perl = TRUE)
fst$chr.name[grepl("Scaffold", fst$chr.name)] <- gsub("(Scaffold_\\d+)__.*", "\\1", fst$chr, perl = TRUE)[grepl("Scaffold", fst$chr.name)]

fst$chr.num <- as.numeric(gsub("(Chr|Scaffold_)(\\d+)(-|_)?.*", "\\2", fst$chr.name, perl = TRUE))
# is it a chromosome or a scaffold
fst$feature <- gsub("(Chr|Scaffold).*", "\\1", fst$chr.name, perl = TRUE)

# make a cumulative position column to concatenate the scaffolds into a single
fst <- fst[order(fst$chr.num, fst$middle.pos), ]

cum.dist <- 0
fst$cum.dist <- 0

for (scf in unique(fst$chr.name)) {
  if (! grepl(pattern = "Scaffold", scf)) {
    # it's a chromosome, keep same positions
    fst$cum.dist[fst$chr.name == scf] <- fst$middle.pos[fst$chr.name == scf]
  } else {
    fst$cum.dist[fst$chr.name == scf] <- cum.dist + fst$middle.pos[fst$chr.name == scf]
    cum.dist <- cum.dist + fst$middle.pos[fst$chr.name == scf][length(fst$middle.pos[fst$chr.name == scf])]
  }
}

# make a new column to plot the scaffold in the same facet
fst$chrFacet <- NA
fst$chrFacet[fst$feature == "Chr"] <- fst$chr.name[fst$feature == "Chr"]
fst$chrFacet[fst$feature == "Scaffold"] <- "Scaffolds"

fst$mb.pos <- fst$cum.dist/1000000

# plot ----------------------------------------------------------------
# whole genome
ggplot(data = fst,
       aes(x = cum.dist/1000000,
           y = value)) +
  geom_point(col = "#737373",
             size = 0.2) +
  geom_smooth(col = "#e6002e") +
  facet_wrap(~chrFacet, 
             nrow = 1,
             scales = "free_x") +
  theme_classic() +
  xlab("Position (Mb)") +
  ylab("Fst value") +
  theme(strip.background = element_blank())
ggsave("figures/exploratory/fst/fst_admixturePops2575_unfolded_wind200kstep100k.png",
       width = 8,
       height = 3,
       dpi = 600)
