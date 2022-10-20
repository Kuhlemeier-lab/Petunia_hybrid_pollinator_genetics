## Custom functions for GWAS analysis
# Marta Binaghi
# 30th September 2020
# Last modified: 09/06/2022

# qq plot
myqqplot <- function(df, p, title) {
  n <- dim(df)[1]
  x <- -log10((1:n)/n)
  y <- rev(sort(-log10(df[ , p])))
  plot(x, y,
       xlab = "Theoretical quantiles",
       ylab = "Observed quantiles",
       main = title)
  abline(0, 1,
         col = "#e6002e")
}

# manhattan plot
manhattan <- function(df, pval, main, threshold = 0, correct) {
  df$chr <- as.integer(gsub(pattern = "^(Chr|Scaffold_)(\\d+)[_-].*",
                            replacement = "\\2",
                            df$rs, 
                            perl = TRUE))
  df$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            df$rs, 
                            perl = TRUE))
  # sort
  df <- df[order(df$chr, df$pos), ]
  # calculate remainder for every chromosome (to assign alternate colours)
  chr_rem <- setNames(1:length(unique(df$chr))%%2,
                      unique(df$chr))
  df$chr_col <- as.character(chr_rem[as.character(df$chr)])
  df$chr_col[df$chr_col == "0"] <- "#d2d2d2"
  df$chr_col[df$chr_col == "1"] <- "#aeaeae"
  # calculate cumulative position
  prev <- 0
  cum.pos <- numeric()
  for (chr in unique(df$chr)) {
    thischr <- df[df$chr == chr, ]
    cum.pos <- c(cum.pos, thischr$pos + prev)
    prev <- cum.pos[length(cum.pos)]
  }
  # save P values in a vector and make the -log10
  p <- -log10(df[ , pval])
  if (correct) {
    # calculate multiple-testing threshold
    thresh_value <- -log10( (threshold) / length(p))
  } else {
    thresh_value <- -log10( threshold )
  }
  # check if threshold fits in plot, otherwise set y axis limit accordingly
  if (thresh_value > max(p)) {
    ytoplim <- thresh_value + ( thresh_value / 10 )
  } else {
    ytoplim <- max(p) + ( max(p) / 10 )
  }
  # plot
  plot(x = cum.pos,
       y = p,
       col = df$chr_col,
       pch = 16,
       ylim = c(0, ytoplim),
       xaxt = "n",
       xlab = "",
       ylab = "-log10(P)",
       main = main)
  if (threshold != 0) {
    abline(h = thresh_value, lty = 2)
  }
}

# manhattan plot for one region
manhattan_region <- function(df, pval, main, threshold_value = 0) {
  df$chr <- as.integer(gsub(pattern = "^(Chr|Scaffold_)(\\d+)[_-].*",
                            replacement = "\\2",
                            df$rs, 
                            perl = T))
  df$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            df$rs, 
                            perl = T))
  # sort
  df <- df[order(df$chr, df$pos), ]
  # save P values in a vector
  p <- df[ , pval]
  # check if threshold fits in plot, otherwise set y axis limit accrdingly
  if (threshold_value > max(p)) {
    ytoplim <- threshold_value + ( threshold_value / 10 )
  } else {
    ytoplim <- max(p) + ( max(p) / 10 )
  }
  # plot
  plot(x = df$pos/1000000,
       y = -log10(p),
       col = "#aeaeae",
       pch = 16,
       ylim = c(0, ytoplim),
       xlab = "Position (Mb)",
       ylab = "-log10(P)",
       main = main)
  if (threshold_value != 0) {
    abline(h = threshold_value, lty = 2)
  }
}




# PIP plot
pip_plot <- function(mydata, title) {
  # add linkage group column (chr)
  chr <- as.numeric(gsub("^(Chr|Scaffold_)(\\d+)[_-].*", "\\2", mydata$rs, perl = T))
  mydata["chr"] <- chr
  # add position column (pos)
  pos <- as.numeric(gsub(".*[-](\\d+)", "\\1", mydata$rs, perl = T))
  mydata["pos"] <- pos
  # sort by linkage group and position
  mydata.sort <- mydata[order(as.numeric(mydata$chr), mydata$pos), ]
  # get list of linkage groups/chromosomes
  chrs <- sort(as.numeric(unique(chr)))
  # calculate remainder for every chromosome (to assign alternate colours)
  chr_rem <- setNames(1:length(chrs)%%2,
                      chrs)
  mydata.sort$chr_col <- as.character(chr_rem[as.character(mydata.sort$chr)])
  mydata.sort$chr_col[mydata.sort$chr_col == "0"] <- "#d2d2d2"
  mydata.sort$chr_col[mydata.sort$chr_col == "1"] <- "#e6002e"
  # calculate cumulative position
  prev <- 0
  cum.pos <- numeric()
  for (chr in chrs) {
    thischr <- mydata.sort[mydata.sort$chr == chr, ]
    cum.pos <- c(cum.pos, thischr$pos + prev)
    prev <- cum.pos[length(cum.pos)]
  }
  # sparse effect size, used for dot size
  z <- mydata.sort$eff
  plot(cum.pos,
       mydata.sort$gamma,
       col = mydata.sort$chr_col,
       pch = 20,
       cex = rescale(x = z, to = c(0.5, 4)),
       xaxt = "n",
       xlab = "",
       ylab = "PIP",
       main = title)
}

# Beta plot
betaplot <- function(df, pval, threshold = 0, correct, main) {
  df$chr <- as.integer(gsub(pattern = "^(Chr|Scaffold_)(\\d+)[_-].*",
                            replacement = "\\2",
                            df$rs, 
                            perl = T))
  df$pos <- as.integer(gsub(pattern = "^.*-(\\d+)$",
                            replacement = "\\1",
                            df$rs, 
                            perl = T))
  # sort
  df <- df[order(df$chr, df$pos), ]
  # calculate remainder for every chromosome (to assign alternate colours)
  chr_rem <- setNames(1:length(unique(df$chr))%%2,
                      unique(df$chr))
  df$chr_col <- as.character(chr_rem[as.character(df$chr)])
  df$chr_col[df$chr_col == "0"] <- "#d2d2d2"
  df$chr_col[df$chr_col == "1"] <- "#e6002e"
  # calculate cumulative position
  prev <- 0
  cum.pos <- numeric()
  for (chr in unique(df$chr)) {
    thischr <- df[df$chr == chr, ]
    cum.pos <- c(cum.pos, thischr$pos + prev)
    prev <- cum.pos[length(cum.pos)]
  }
  # save P values in a vector
  p <- df[ , pval]
  if (correct) {
    # calculate multiple-testing p value
    p <- p*length(p)
  }
  df$cum.pos <- cum.pos
  df$p <- p
  # subset to keep only sites passing p-value filter
  if (threshold != 0) {
    df_plot <- df[df$p <= threshold, ]
  } else {
    df_plot <- df
  }
  # plot
  plot(x = df_plot$cum.pos,
       y = df_plot$beta,
       col = df_plot$chr_col,
       pch = 16,
       xaxt = "n",
       xlab = "",
       ylab = "Beta",
       main = main)
  # add errorbars
  arrows(df_plot$beta - df_plot$se, 
         df_plot$cum.pos, 
         df_plot$beta + df_plot$se, 
         df_plot$cum.pos, 
         angle = 90, 
         length = 0.05, 
         code = 3)
}


# fst plot
plotFst <- function(df, main, threshold = 0) {
  df$chrname <- as.integer(gsub(pattern = "^(Chr|Scaffold_)(\\d+)[_-].*",
                            replacement = "\\2",
                            df$chr, 
                            perl = T))
  # sort
  df <- df[order(df$chrname, df$middle.pos), ]
  # calculate remainder for every chromosome (to assign alternate colours)
  chr_rem <- setNames(1:length(unique(df$chrname))%%2,
                      unique(df$chrname))
  df$chr_col <- as.character(chr_rem[as.character(df$chrname)])
  df$chr_col[df$chr_col == "0"] <- "#d2d2d2"
  df$chr_col[df$chr_col == "1"] <- "#e6002e"
  # calculate cumulative position
  prev <- 0
  cum.pos <- numeric()
  for (chr in unique(df$chrname)) {
    thischr <- df[df$chrname == chr, ]
    cum.pos <- c(cum.pos, thischr$middle.pos + prev)
    prev <- cum.pos[length(cum.pos)]
  }
  # calculate threshold (ie quantile)
  thresh_value <- quantile(df$value, probs = threshold, na.rm = TRUE)
  # check if threshold fits in plot, otherwise set y axis limit accordingly
  if (thresh_value > max(df$value, na.rm = TRUE)) {
    ytoplim <- thresh_value + ( thresh_value / 10 )
  } else {
    ytoplim <- max(df$value, na.rm = TRUE) + ( max(df$value, na.rm = TRUE) / 10 )
  }
  # plot
  plot(x = cum.pos,
       y = df$value,
       col = df$chr_col,
       pch = 16,
       ylim = c(0, ytoplim),
       xaxt = "n",
       xlab = "",
       ylab = "Fst",
       main = main)
  if (threshold != 0) {
    abline(h = thresh_value, lty = 2)
  }
}
