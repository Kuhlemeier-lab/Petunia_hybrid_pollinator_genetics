# Admixture analysis, natural hybrids
# Marta Binaghi
# Last modified 11/05/2022
# 
# The admixture analysis was performed with ngsadmix, the scripts used were:
# 11a_ngsadmix.sh to run the admixture analysis
# The dataset used was therefore: hardfiltered_biallelic_cr09_mm005.beagle.gz.
# Note that in this file the samples are numbered from 0 to 69. The correspondence between
# the beagle sample (0-69) and the raw reads ID is shown in file data/reads_sample_ID.csv.
# I tested K from 1 to 8, repeating each K 10 times. I then select one run per 
# K (the one with the best likelihood) and apply the Evanno's method to define
# the most likely K.

# libraries ---------------------------------------------------------------
library(data.table)     # allows subsetting df based on values in another df
library(dplyr)
library(Hmisc)          # add error bars on plot
library(RColorBrewer)

wdir <- "hybrids"
setwd(wdir)


# Find best run per each K ------------------------------------------------
# I have parsed the log files of ngsadmix using a custom python script, 
# ngsadmix_outparser.py
# This produces the log_stats.csv file, which I read into R. 
# The best likelihood in each K tested is found.

logstat <- read.table("data/raw/admixture/log_stats.csv",
                      sep = ",",
                      dec = ".",
                      stringsAsFactors = F,
                      header = T)

# remove NAs
logstat <- na.omit(logstat)
# sort by K and rep
logstat <- logstat[order(logstat$K, logstat$rep), ]

# find run with highest likelihood per K
best_reps <- numeric()
ks <- numeric()
for ( k in levels(as.factor(logstat$K)) ) {
  ks <- c(ks, as.numeric(k))
  kdf <- logstat[logstat$K == k, ]
  best_rep <- kdf$rep[kdf$likelihood == max(kdf$likelihood)][1]
  best_reps <- c(best_reps, best_rep)
}
best_reps_df <- as.data.frame(cbind(ks, best_reps))
colnames(best_reps_df) <- c("K", "rep")

keys <- c("K", "rep")
tData <- data.table(logstat, key=keys)
tBest <- data.table(best_reps_df, key=keys)

best <- tData[tBest, ]
write.csv(best,
          file = "data/raw/admixture/best_likelihood_run.csv",
          row.names = F,
          quote = F)


# most likely K -----------------------------------------------------------
# I apply Evanno's method to find the most likely K.
ksummary <- logstat %>%
  group_by(K) %>%
  summarise(meanL = mean(likelihood),
            sdL = sd(likelihood)
  )
           

deltak <- numeric()
for (i in 2:10) {
  delta <- abs(ksummary$meanL[ksummary$K == (i+1)] -
        (2 * ksummary$meanL[ksummary$K == i] ) +
        ksummary$meanL[ksummary$K == (i-1)]
        ) / ksummary$sdL[ksummary$K == i]
  deltak <- c(deltak, delta)
}

# And plot the log likelihood
plot(x = ksummary$K,
     y = ksummary$meanL,
     type = "n")
pdf(file = "figures/exploratory/admixture/plot_evanno.pdf",
    width = 6.5,
    height = 4)
par(mfrow = c(1, 2))

with(data = ksummary,
     expr = errbar(K, meanL, meanL+sdL, meanL-sdL, add=F, pch=1, cap=.015,
                   ylab = "Mean log-likelihood")
)
title(main = "Log-likelihood mean and SD")
# plot delta K
plot(x = 2:7,
     y = deltak,
     main = "DeltaK",
     xlab = "K",
     ylab = "deltaK")
dev.off()
# The most likely K is 2.

# Admixture proportion plotting -------------------------------------------
# Now that the K has been defined, I can plot the proportion of admixture per
# each sample. To do this I wrote a couple of functions.

# admixture plot functions -------------------------------------------------
read_meanQ <- function(file) {
  qtab <- read.table(file,
             sep = "",
             dec = ".",
             header = F,
             stringsAsFactors = F)
  Ks <- dim(qtab)[2]
  colnames(qtab) <- paste0("Q", Ks, ".", 1:Ks)
  return(qtab)
}
order_dataOLD <- function(qmat) {
  # check if label column does not exists
  if (!colnames(qmat)[1] %in% c("label", "labels")) {
    # assign individuals to populations
    pop_ass <- apply(qmat, 1, which.max)
    qmat <- cbind(label = pop_ass,
                  qmat[ ])
  }
  # order data
  qmat <- qmat[do.call(order, qmat), ]
  return(qmat)
}
order_data <- function(qmat) {
  # round admixture values so that the sorting works correctly
  nonsorted_ids <- 1:nrow(qmat)
  rounded <- qmat
  rounded[ , sapply(rounded, is.numeric)] <- apply(rounded[ , sapply(rounded, is.numeric)],
                                                   c(1, 2),
                                                   round, 
                                                   digits = 2)
  # check if label column does not exists
  if (!colnames(rounded)[1] %in% c("label", "labels")) {
    # assign individuals to populations
    pop_ass <- apply(rounded, 1, which.max)
    rounded <- cbind(label = pop_ass,
                     rounded[ ])
  }
  # order data
  rounded$sorting_ids <- nonsorted_ids
  rounded <- rounded[do.call(order, rounded), ]
  # use the sorted and rounded q matrix to sort the original (non rounded) matrix
  qmat <- qmat[match(rounded$sorting_ids, nonsorted_ids), ]
  return(qmat)
}
structure_plot <- function(mydata) {
  lab_exist <- F
  # get sample id
  id <- mydata$id
  # check if first column is labels
  if (colnames(mydata)[1] %in% c("label", "labels")) {
    labels <- as.factor(mydata[ , 1])
    # drop unnecessary cols from admixture data
    admix_data <- subset(mydata, select = -c(label, id))
    #admix_data <- mydata[ , -1]
    # logical to know if needs plotting of label groups
    lab_exist <- T
  } else {
    admix_data <- subset(mydata, select = -c(id))
  }
  # calculate K
  K <- dim(admix_data)[2]
  # plot
  #print(head(t(as.matrix(admix_data))))
  #print(admix_data$id)
  myplot <- barplot(t(as.matrix(admix_data)),
                    col = brewer.pal(n = K, name = "Paired"),
                    cex.names = .7,
                    las = 2,
                    ylim = c(0, 1.3),
                    yaxt = "n",
                    xlab = "ID of individuals",
                    ylab = "Ancestry proportion",
                    names.arg = id)
  myplot
  axis(2,
       at = seq(0, 1, 0.2),
       labels = seq(0, 1, 0.2),
       las = 2)
  if (lab_exist) {
    # calculate start and end of groups, to draw lines
    labels <- as.factor(mydata$label)
    line_starts <- numeric()
    line_ends <- numeric()
    for (label in levels(labels)) {
      start <- myplot[which(labels == label)[1]]
      end <- myplot[which(labels == label)[length(which(labels == label))]]
      line_starts <- c(line_starts, start)
      line_ends <- c(line_ends, end)
    }
    # add custom top "axis"
    segments(x0 = line_starts - 0.2,
             y0 = 1.05,
             x1 = line_ends + 0.2,
             y1 = 1.05,
             col = brewer.pal(n = length(line_starts), name = "Dark2"))
    text(x = ((line_ends - line_starts) / 2 ) + line_starts,
         y = 1.2,
         labels = levels(labels),
         srt = 90)
  }
}

# I import the Q matrices.
myQmatrix <- "data/raw/admixture/hf_ba_cr09_mm005_k"
for (K in seq(2, 8) ) {
  dataset_name <- paste0("qmat", K)
  dataset_file <- paste0(myQmatrix,
                         K,
                         "_run",
                         best$rep[best$K == K],
                         ".qopt")
  assign(dataset_name,
         read_meanQ(dataset_file)
  )
}

# I then import the labels and IDs
# import custom labels
labels <- read.table("data/clean/reads_sample_ID.csv",
                     sep = ",",
                     header = TRUE,
                     stringsAsFactors = FALSE)
# keep only one lane in the metadata because we only need the sample information
labels <- labels[labels$lane == 1, ]
labels$batch <- 2018
labels$batch[grepl("KMH", labels$sra_library_id)] <- 2016
labels <- labels[order(labels$id_gl), ]

# Then I order the data and make one pdf per K {2-8}.
# In each pdf there are three plots, showing the same admixture values, but with
# different order of the individuals, to visualise different groups: by location (Ponto 143 and Pedra
# da Cruz), by parental plant, by admixture group.
# The pdfs are stored in figures/exploratory/admixture.

# data sorting is done first on the label column, then on the rest.
# if labels are not present, I assign every individual to a population by
# checking its bigger admixture component.
for (K in seq(2, 8) ) {
  pdf(file = paste0("figures/exploratory/admixture/K", K, ".pdf"),
      title = "Ancestry proportions of hybrid populations",
      width = 6.5 * 1.3,
      height = 8.66 * 1.3)
  par(mfrow=c(3,1))
  dataset_name <- paste0("qmat", K)
  # make three plots, ordering by inferred population, parental family,
  #  and location (Pedra da Cruz or Ponto 143)
  label_column <- c("location", "parental_id")
  for (i in 1:length(label_column)) {
    to_order <- cbind(label = labels[[label_column[i]]],
                      get(dataset_name),
                      id = labels$new_id_short)
    to_plot <- order_data(to_order)
    structure_plot(to_plot)
    title(main = paste0("Ancestry proportions for K = ",
                        K,
                        ", groups defined by ",
                        label_column[i]))
  }
  to_plot <- order_data(get(dataset_name))
  structure_plot(to_plot)
  title(main = paste0("Ancestry proportions for K = ",
                      K,
                      ", groups defined by major ancestry proportion"))
  dev.off()
}


# export admixture matrices -----------------------------------------------
# I make a dataframe containing all the metadata of the samples (ID, parental, batch, ...)
# and add the admixture values per each K tested, per each sample.
# save to which K group each ID belongs in a dataframe
full_df <- cbind(labels,
      qmat2,
      qmat3,
      qmat4,
      qmat5,
      qmat6,
      qmat7,
      qmat8,
      stringsAsFactors = F)

write.csv(full_df,
          file = "data/raw/admixture/admixture_matrices.csv",
          quote = F,
          row.names = F)

# I also save K2 in a separated file.
# save only K2 groups
q2_group <- rep(2, times = 70)
q2_group[qmat2$Q2.1 >= 0.50000000] <- 1
q2 <- data.frame(id = labels$new_id_short,
                 id_gl = labels$id_gl,
                 location = labels$location,
                 q2_group = q2_group,
                 Q2.1 = qmat2$Q2.1,
                 Q2.2 = qmat2$Q2.2,
                 stringsAsFactors = FALSE)
write.csv(q2,
          file = "data/raw/admixture/q2_admixture.csv",
          quote = F,
          row.names = F)

# save the IDs of individuals in group <= 0.25 and >=0.75 for later use in the
# divergence analysis
# admixture group 1, adm >= 0.75
admixture1 <- "admixture,1,"
admixture1members <- paste0(q2$id[q2$Q2.1 >= 0.75], collapse = ",")
admixture1str <- paste0(admixture1, admixture1members, collapse = "")
# admixture group 2, adm <= 0.25
admixture2 <- "admixture,2,"
admixture2members <- paste0(q2$id[q2$Q2.1 <= 0.25], collapse = ",")
admixture2str <- paste0(admixture2, admixture2members, collapse = "")

myfile <- "data/raw/samples_in_admixture_groups.csv"
cat(admixture1str, file = myfile, append = FALSE, sep = "\n")
cat(admixture2str, file = myfile, append = TRUE, sep = "\n")


# manuscript plot for K = 2 -----------------------------------------------
# And the final plot for K = 2, with red colour showing admixture group 2, and white group 1.
q2 <- cbind(qmat2,
            id = labels$new_id_short)
mymat <- order_data(q2)
# for consistency among the figures, I invert the order of the q2 matrix,
# so that first comes group 1 (ie Q2.1 >= 0.50) and then group 2 (ie Q2.1 < 0.50)
mymat <- mymat[rev(rownames(mymat)), ]
mypal <- c('#e6002e', '#d9d9d9')

pdf(file = "figures/exploratory/admixture/admixture_k2.pdf",
    width = 10,
    height = 4,
    title = "Admixture proportions, K2")
myplot <- barplot(t(as.matrix(mymat[ , 2:3])),
                  col = mypal,
                  las = 2,
                  cex.names = .7,
                  ylim = c(0, 1),
                  xlab = "ID of individuals",
                  ylab = "Ancestry proportion",
                  names.arg = mymat$id)
dev.off()

#Note that for K = 2, individuals group as follow:
sum(qmat2$Q2.1 >= 0.75) # individuals with >= 0.75 group 1
sum(qmat2$Q2.1 <= 0.25) # individuals with <= 0.25 group 2
sum(qmat2$Q2.1 < 0.75 & qmat2$Q2.1 > 0.25) # individuals with admixture between 0.25 and 0.75
sum(qmat2$Q2.1 >= 0.50) # individuals with >= 0.50 group 1
sum(qmat2$Q2.1 < 0.50) # individuals with < 0.50 group 2
