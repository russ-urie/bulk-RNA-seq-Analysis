# File Info ----
# Author: Russell Ricks Urie, PhD
# Title: Skin Transplant Training Cohort RNAseq Analysis
# Date Created: Jan 2021    Date Last Updated: 03-10-2022
# Info: Bulk RNA-seq of Skin Transplant Training Cohort from Aaron/Diana
sample_set <- "STx_TC"
# 8 skin transplants (STx) performed, 3 good allogeneic & syngeneic. 3 scaffolds 
# were implanted and removed at day 0, day 7, and day 13 after T cell transfer.
currentDate <- Sys.Date()
# Sources Link: https://docs.google.com/document/d/1QVrDYwP-M40FaVOiZ3BCLiY9tCoB8Ra9ZLaMU-g6mYA/edit?usp=sharing)

# Load libraries ----
# install.packages("pacman") # Be sure to install pacman if you haven't before
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, glmnet,
  colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, apeglm, rgl, qpcR, boot, caret, ggvenn,
  devtools, reshape2, gridExtra, factoextra, edgeR, cowplot, pheatmap, coefplot, grid,
  randomForest, ROCR, genefilter, Hmisc, rdist, factoextra, ggforce, NormqPCR, ggpubr)

# A) PREPARE DATA ----------------------------------------------------------
  # Raw and CPM data are read in & organized. Raw data is pseudo-normalized.
# 1 Organize Data ----
samples <- c("Sd00m1","Ad07m5","Ad07m6","Ad07m8","Sd13m1","Sd13m3","Sd13m4","Ad13m5","Ad13m6",
             "Ad13m8","Sd00m3","Sd00m4","Ad00m5","Ad00m6","Ad00m8","Sd07m1","Sd07m3","Sd07m4")
sample_num <- length(samples)
re_order <- c(1,11,12,16,17,18,5,6,7,13,14,15,2,3,4,8,9,10) # Proper order
samples <- samples[re_order] # Samples put in order
r_file <- read.table(".\\Inputs\\STX_TC\\counts_raw.txt", sep = "", header=T)
r_file <- r_file[ grep("Rik", r_file$X, invert = T) , ] # Remove Riken genes
genes <- make.names(r_file[,1], unique = T) # Cant start w/ numbers as row names
r_counts <- r_file[,-1] # Remove first column of gene names temporarily
r_counts <- r_counts[,re_order] # Counts put in order
colnames(r_counts) <- c(samples) # Add sample names as column names
rownames(r_counts) <- genes # Add gene names as row names
r_counts <- r_counts[rowSums(r_counts)>0,] # Remove 0 expression genes
genes <- rownames(r_counts) # Non-zero genes
n_file <- read.table(".\\Inputs\\STX_TC\\counts_per_million.txt", sep="", header=T)
n_file <- n_file[ grep("Rik", n_file$X, invert = T) , ] # Remove Riken genes
n_counts <- n_file[,-1] # Remove first column of gene names temporarily
n_counts <- n_counts[,re_order] # Samples now back in order!
colnames(n_counts) <- c(samples) # Add sample names as column names
n_counts <- n_counts[rowSums(n_counts)>0,] # Remove 0 expression genes
rownames(n_counts) <- genes # Add gene names back as row names

# 2 Define Factors ----
cohort <- c("SynDay00","SynDay07","SynDay13","AlloDay00","AlloDay07","AlloDay13")
donor <- c(rep("Syngeneic", 9), rep("Allogeneic", 9))
time <- c("Day00", "Day07", "Day13") # As strings for labeling plots
day <- c(0,7,13) # As numeric
mouse <- c(rep(c("M1", "M3", "M4"), 3), rep(c("M5", "M6", "M8"),3))
donor_num <- 2 # Used to shorten later code
time_num <- 3
mouse_num <- 3
cohort <- rep(cohort, times=1, each=mouse_num) # Used to categorize samples
time <- rep(time, times=donor_num, each=mouse_num)
immune <- c(rep("Healthy",3), rep("Immunized",6), rep("Healthy",3), rep("Immunized",6))
day <- rep(day, times=donor_num, each=mouse_num)
grouped <- c(rep("Healthy", 12), rep("Early_Rejecting", 3), rep("Late_Rejecting", 3))
simple <- c(rep("Healthy", 12), rep("Rejecting", 6))
clin_info <- cbind.data.frame(donor, cohort, time, mouse, immune, simple, grouped)
rownames(clin_info) <- samples
clin_info <- cbind(clin_info, c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1)) # 0=healthy
colnames(clin_info)[8] <- "state"

# 3 Pseudo-normalize ----
ps_counts <- log(n_counts + 1, base = 2)
ps_LF <- melt(ps_counts, variable.name="Samples", value.name="Count") # Long form
# Use first 4 letters of sample names as conditions for labeling plots
ps_LF <- data.frame(ps_LF, Condition = substr(ps_LF$Samples, 1, 4))
ps_LF$Condition <- factor(ps_LF$Condition, # Order groups, not alphabetical
                            levels=c("Sd00", "Sd07", "Sd13", "Ad00", "Ad07", "Ad13"))
# Remove genes when 3/4 of samples have 0 expression
ps_counts <- ps_counts[rowSums(ps_counts == 0) <= sample_num*(3/4),]
genes <- rownames(ps_counts)
r_counts <- r_counts[genes,]
n_counts <- n_counts[genes,]

# 4 Data Overview ----
# To assess the data quality; can be skipped after first pass.
# Raw count & Pseudo-normalized histograms
if (sample_num <= 20) { # If 20 samples or less, print to 1 page
  par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
} else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
par(mar = c(4,1,2,1), oma = c(0.5,1,0,0)) # Adjust margins to better use space
sapply(1:sample_num, function(x) hist(n_counts[,x], ylim=c(0,1000), breaks=100, 
                                      xlab=samples[x], ylab = "", main = ""))
mtext("Counts per Million, 100 bins", side = 3, line = -2, outer = T) # Add title
if (sample_num <= 20) { # If 20 samples or less, print to 1 page
  par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
} else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
par(mar = c(4,1,2,1), oma = c(0.5,1,0,0)) # Adjust margins to better use space
sapply(1:sample_num, function(x) hist(ps_counts[,x], ylim=c(0,1000), breaks=100, 
                                      xlab=samples[x], ylab = "", main = ""))
mtext("log[2] (CPM + 1), 100 bins", side = 3, line = -2, outer = T) # Add title
# Box & Density Plots, how samples compare in mean, variance, & potential outliers
p_pseudobox <- ggplot(ps_LF, aes(x=Samples, y=Count, fill=Condition)) + 
  geom_boxplot() + xlab("") + ylab(expression(log[2](count+1))) +
  scale_x_discrete(breaks=samples, labels=mouse) + ggtitle("Pseudo Unfilt") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
p_pseudobox
p_pseudodens <- ggplot(ps_LF, aes(x=Count, color=Samples, fill=Samples)) + 
  geom_density(alpha=0.2, size=0.75) + facet_wrap(~ Condition) + 
  theme(legend.position="top") + xlab(expression(log[2](count + 1))) + 
  ggtitle("Pseudo Unfilt")
p_pseudodens
p_samplediff <- sample.diff(ps_counts, samples) # Likeness of samples
# p_MAplots <- plot.topMA(ps_counts) # Sample likeness in detail (VERY SLOW)
ggbiplot(prcomp(t(ps_counts), scale.=T), ellipse=T, groups=time, var.axes=F, 
         labels=mouse, var.scale=1, circle=T) + ggtitle("Pseudo Unfilt") + 
  theme_classic() 
ggscreeplot(prcomp(t(ps_counts), scale.=T))
par(mfrow=c(1,1))

# GSEA Unfilt ----
# fac_DEseq <- grouped
# coldata <- data.frame(cbind(fac_DEseq))
# row.names(coldata) <- samples # DESeq2 uses raw counts; rows:genes & cols:samples
# dds <- DESeqDataSetFromMatrix(countData=r_counts, colData=coldata, design = ~ fac_DEseq)
# paste(nrow(dds), " genes input into DESeq2 for GSEA", sep="")
# dds <- DESeq(dds)
# human_genes1 <- convertMouseGeneList(as.array(rownames(get_counts(dds))))
# human_genes1 <- human_genes1[!duplicated(human_genes1[1]),]
# human_genes2 <- as.data.frame(human_genes1[,2])
# row.names(human_genes2) <- human_genes1[,1]
# colnames(human_genes2) <- "human_genes"
# combined_list <- merge(get_counts(dds), human_genes2, by = "row.names", all.x = T)
# nulllist <- combined_list[is.na(combined_list$human_genes),]
# nulllist$human_genes <- nulllist$Row.names
# notnulllist <- combined_list[!is.na(combined_list$human_genes),]
# combined_list2 <- rbind(nulllist, notnulllist)
# gct_format <- cbind(combined_list2[,length(combined_list2)], NA, 
#                     combined_list2[,3:length(combined_list2)-1])
# gct_format <- rbind("", "", colnames(gct_format), gct_format)
# gct_format[1,1] <- "#1.2"
# gct_format[2,1] <- nrow(combined_list2)
# gct_format[2,2] <- ncol(combined_list2)-2
# gct_format[3,1] <- "NAME"
# gct_format[3,2] <- "Description"
# write.table(gct_format, file = paste(currentDate, sample_set, "_unfilt.gct", sep=""), 
#             quote=F, row.names=F, col.names=F, sep ="\t")
# cls_format <- rbind("", "", colnames(combined_list2[3:ncol(combined_list2)-1]))
# cls_format[3,] <- fac_DEseq
# groups1 <- fac_DEseq[!duplicated(fac_DEseq)]
# cls_format[1,1:3] <- c(ncol(cls_format), length(groups1),1)
# cls_format[2,1:(length(groups1)+1)] <- c("#", groups1)
# write.table(cls_format, file = paste(currentDate,sample_set, "_unfilt_classes.cls", 
#                                      sep=""), quote=F, row.names=F, col.names=F)

# B) REDUCING NOISE --------------------------------------------------------
# 5 Variance Cutoffs ----
# Genes with at least moderate variance among all samples but low healthy variance
par(mfrow=c(1,1)) # Return to default page orientation
par(mar=c(4,4,2,1))
mean_cut <- 0.2
var_cut <- 0.01
df <- cbind.data.frame("x"=rowMeans(ps_counts), "y"=rowVars(ps_counts)/rowMeans(ps_counts))
g_cutoff <- rownames(subset(df, x > mean_cut & y > var_cut))
with(df, plot(x, y, pch=20, cex=1, main=paste("All Var vs Mean, ", length(g_cutoff), 
              " genes w/ Mean > ", mean_cut, " & ", "Var > ", var_cut, sep=""), 
              xlab="Mean across 18 samples", ylab="Var/Mean across 18 samples"))
with(subset(df, x > mean_cut & y > var_cut), points(x, y, pch=20, col="red3", cex=1))
abline(v=mean_cut, col="red3", lty=2, lwd=1.5) # Line for Mean cutoff
abline(h=var_cut, col="red3", lty=2, lwd=1.5) # Line for Var cutoff
ggbiplot(prcomp(t(ps_counts[ g_cutoff,]), scale.=T), ellipse=T, groups=simple, 
         var.axes=F, labels=mouse, var.scale=1, circle=T) + theme_classic() +
  ggtitle("Pseudo: High Var")
ps_h <- ps_counts[ , simple == "Healthy"]
mean_cut <- 0.1
var_cut <- 0.8
df <- cbind.data.frame("x"=rowMeans(ps_h), "y"=rowVars(ps_h)/rowMeans(ps_h))
g_h <- rownames(subset(df, x >= mean_cut & y <= var_cut))
with(df, plot(x, y, pch=20, cex=1, main=paste("Healthy Var vs Mean, ", length(g_h), 
              " genes with Mean > ", mean_cut, " & ", "Var > ", var_cut, sep=""), 
              xlab="Mean across 12 H samples", ylab="Var/Mean across 12 H samples"))
with(subset(df, x >= mean_cut & y <= var_cut), points(x, y, pch=20, col="red3", cex=1))
abline(v=mean_cut, col="red3", lty=2, lwd=1.5) # Line for Mean cutoff
abline(h=var_cut, col="red3", lty=2, lwd=1.5) # Line for Var cutoff
g1 <- Reduce(intersect, list(g_cutoff, g_h))
ggbiplot(prcomp(t(ps_counts[ g1,]), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=mouse, var.scale=1, circle=T) + theme_classic() +
  ggtitle("Pseudo: High Var, Low Healthy Var")
ps_r <- ps_counts[ , simple == "Rejecting"]
mean_cut <- 0.1
var_cut <- 0.8
df <- cbind.data.frame("x"=rowMeans(ps_r), "y"=rowVars(ps_r)/rowMeans(ps_r))
g_r <- rownames(subset(df, x >= mean_cut & y <= var_cut))
with(df, plot(x, y, pch=20, cex=1, main=paste("Rej Var vs Mean, ", length(g_r), 
              " genes with Mean > ", mean_cut, " & ", "Var > ", var_cut, sep=""), 
              xlab="Mean across 6 Rej samples", ylab="Var/Mean across 6 Rej samples"))
with(subset(df, x >= mean_cut & y <= var_cut), points(x, y, pch=20, col="red3", cex=1))
abline(v=mean_cut, col="red3", lty=2, lwd=1.5) # Line for Mean cutoff
abline(h=var_cut, col="red3", lty=2, lwd=1.5) # Line for Var cutoff
g1 <- Reduce(intersect, list(g_cutoff, g_h, g_r))
ggbiplot(prcomp(t(ps_counts[ g1,]), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=mouse, var.scale=1, circle=T) +   theme_classic() +
  ggtitle("Pseudo: High Var, Low Healthy Var, Low Rejecting Var")
ps_h <- ps_counts[ , simple == "Healthy"]
ps_r <- ps_counts[ , simple == "Rejecting"]
mean_cut <- 0.7
var_cut <- 0.7
df <- cbind.data.frame("x"=rowVars(ps_h)/rowMeans(ps_h), 
                       "y"=rowVars(ps_r)/rowMeans(ps_r))
g_r <- rownames(subset(df, x <= mean_cut & y <= var_cut))
with(df, plot(x, y, pch=20, cex=1, main=paste("Grouped Var vs Mean, ", length(g_r), 
    " genes with Var > ", var_cut, sep=""), 
    xlab="Var/Mean across 12 Healthy samples", ylab="Var/Mean across 6 Rej samples"))
with(subset(df, x <= mean_cut & y <= var_cut), points(x, y, pch=20, col="red3", cex=1))
abline(v=mean_cut, col="red3", lty=2, lwd=1.5) # Line for Mean cutoff
abline(h=var_cut, col="red3", lty=2, lwd=1.5) # Line for Var cutoff
g1 <- Reduce(intersect, list(g_cutoff, g_r))
ggbiplot(prcomp(t(ps_counts[ g_r,]), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=mouse, var.scale=1, circle=T) + theme_classic() +
  ggtitle("Pseudo: High Var, Low Healthy Var, Low Rejecting Var")
ps_var <- ps_counts[ g1, ]
n_var <- n_counts[ g1, ]
r_var <- r_counts[ g1, ]
ps_counts <- ps_var

# 6 Assess Many Filters ----
# Basically the same thing as running DESeq2, just stripped down.
# filt_low <- seq(1,6.5,0.5)
# filt_high <- seq(16,12,-2)
# g_filt <- filters.genes(ps_counts, filt_low, filt_high)
# g_DE <- filters.DEGs(g_filt, simple)
# num_DEGs <- cbind.data.frame("DEGs"=sapply(1:length(g_DE$DEGs), function(x) 
#   length(g_DE$DEGs[[x]])))
# num_DEGs <- cbind.data.frame("DEGs"=num_DEGs, "Filters"=names(g_DE$DEGs))
# coul <- colorRampPalette(brewer.pal(4, "Greens"))(nrow(num_DEGs))
# p_DEGfilts <- ggplot(data=num_DEGs, aes(x=Filters, y=DEGs)) + 
#   geom_bar(stat="identity", color=coul, fill=coul) + 
#   geom_text(aes(label=DEGs), position=position_dodge(width=0.9), vjust=0) +
#   coord_flip() + ggtitle("DEGs at various filters")
# p_DEGfilts
# filters_SVD <- filters.SVDs(g_filt, simple)
# filters_RF <- filters.RF(g_filt, simple, clin_info, "Rejecting")
# filters_scores <- filters.scores(filters_SVD, filters_RF)
# filters_scoreplots <- filters.scoreplots(filters_scores, grouped, donor, time)
# p_scores <- ggarrange(plotlist=filters_scoreplots$ScorePlots, common.legend=T, 
#                       legend="right", nrow=2, ncol=2)
# p_scores

# 7 Expression Filters ----
coul <- c(rep("red",3), rep("blue",3), rep("green",3), rep("brown",3), 
          rep("orange",3), rep("purple",3)) # Colors for lines
types <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3) # Solid or dashed lines
num_breaks <- 400
lst <- lapply(1:sample_num, function (x) { # Store list of histogram densities
  res <- hist(ps_counts[,x], breaks=num_breaks, plot=F)
  res <- cbind(res$breaks[1:num_breaks], res$counts[1:num_breaks])})
plot.variation(lst, types, samples) # Plot all of the histogram density lines
p_denslines <- recordPlot()
plot.error(lst,types) # Sum of variation from the mean
p_sumerror <- recordPlot()
par(mfrow=c(1,1)) # Return to default page orientation
# This section cannot be skipped!
filt_low <- 2 # These three filters are per sample per gene
filt_group <- 0
filt_high <- 10
keep1 <- rowSums(ps_counts) > (filt_low*sample_num) # Low count filter, required
ps_filt <- ps_counts[keep1,] # Apply filter to pseudonorm counts
dim(ps_filt) # Dimensions after applying low filter
keep2 <- rowSums(ps_filt >= (filt_group)) >= mouse_num # Group filter, optional
ps_filt <- ps_filt[keep2,] # Apply filter to pseudonorm counts
dim(ps_filt) # Dimensions after applying group filter
keep3 <- rowSums(ps_filt) < (filt_high*sample_num) # High count filter, optional
ps_filt <- ps_filt[keep3,] # Apply filter to pseudonorm counts
dim(ps_filt) # Dimensions after applying high filter
ps_high <- ps_counts[rowSums(ps_counts) >= (filt_high*sample_num),]
genes_high <- rownames(ps_high) # To check high filter genes for relevance
g_filt <- rownames(ps_filt) # Vector of filtered genes
r_filt <- r_counts[g_filt,] # Used for DESeq2, VERY IMPORTANT
n_filt <- n_counts[g_filt,] # Used for plots and calcs, important

# 8 Filt Data Overview ----
# Histograms
if (sample_num <= 20) { # If 20 samples or less, print to 1 page
  par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
} else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
par(mar = c(4,1,2,1), oma = c(0.5,1,1,0)) # Adjust margins to better use space
lapply(1:sample_num, function(x) {
  old <- lapply(1:sample_num, function(x) hist(ps_counts[,x], breaks=100, plot=F))
  new <- lapply(1:sample_num, function(x) hist(ps_filt[,x], breaks=100, plot=F))
  plot(old[[x]], ylim=c(0,1000), xlab=samples[x], ylab = "", main = "", 
       border = "gray", col = "gray")
  plot(new[[x]], add=T, border = "green4", col = "green4")})
mtext(paste("Pseudo Filtered:", filt_low, "&", filt_high, 
            "filters"), side = 3, line = -1, outer = T) # Add title
p_pseudofhist <- recordPlot()
par(mfrow=c(1,1)) # Return to default page orientation
par(mar = c(4,4,4,4), oma = c(1,1,1,1))
#  Box Plot, convert to "long form", then plot as in section 4
ps_fLF <- melt(ps_filt, variable.name="Samples", value.name="Count")
ps_fLF <- data.frame(ps_fLF, Condition = substr(ps_fLF$Samples, 1, 4))
ps_fLF$Condition <- factor(ps_fLF$Condition, levels = c("Sd00", "Sd07", "Sd13", 
                                                        "Ad00", "Ad07", "Ad13"))
p_pseudofbox <- ggplot(ps_fLF, aes(x=Samples, y=Count, fill=Condition)) + 
  geom_boxplot() + xlab("") + ylab("") +
  scale_x_discrete(breaks=samples, labels=mouse) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + ylim(0,16) + 
  ggtitle(paste("Pseudo Filtered:", filt_low, "&", filt_high, "filters"))
p_pseudofbox <- ggarrange(p_pseudobox, p_pseudofbox, common.legend=T, 
                          legend="right", nrow=1, ncol=2)
p_pseudofbox
p_pseudofdens <- ggplot(ps_fLF, aes(x=Count, color=Samples, fill=Samples)) + 
  geom_density(alpha=0.2, size=0.75) + facet_wrap(~ Condition) + 
  theme(legend.position="top") + xlab(expression(log[2](count + 1))) + 
  ggtitle(paste("Pseudo Filtered:", filt_low, "&", filt_high, "filters"))
p_pseudofdens <- ggarrange(p_pseudodens, p_pseudofdens, common.legend=T, 
                           legend="top", nrow=1, ncol=2)
p_pseudofdens
p_samplefdiff <-  sample.diff(ps_filt, samples)
p_fPCA <- ggbiplot(prcomp(t(ps_filt), scale.=T), ellipse=T, groups=cohort, 
                   var.axes=F, labels=mouse, var.scale=1, circle=T) + 
  ggtitle("Pseudo Filt") + theme_classic() 
p_fPCA
ggscreeplot(prcomp(t(ps_filt), scale.=T))
# Histogram Densities
coul <- c(rep("red",3), rep("blue",3), rep("green",3), rep("brown",3), 
          rep("orange",3), rep("purple",3)) # Colors for lines
types <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3) # Solid or dashed lines
num_breaks <- 400
lst <- lapply(1:sample_num, function (x) { # Store list of histogram densities
  res <- hist(ps_filt[,x], breaks=num_breaks, plot=F)
  res <- cbind(res$breaks[1:num_breaks], res$counts[1:num_breaks])})
plot.variation(lst, types, samples, ylim1=500, ylim2=500) 
p_denslines <- recordPlot()
plot.error(lst,types) # Sum of variation from the mean
p_sumerror <- recordPlot()
par(mfrow=c(1,1))
if (sample_num <= 20) { # If 20 samples or less, print to 1 page
  par(mfrow=c(4,ceiling(sample_num/4))) # All histograms on a page w/ 4 rows
} else {par(mfrow=c(3,5))} # Histograms on multiple pages w/ 3 rows & 5 cols
par(mar = c(4,1,2,1), oma = c(0.5,1,1,0)) # Adjust margins to better use space
lapply(1:sample_num, function(x) {
  old <- lapply(1:sample_num, function(x) hist(ps_counts[,x], breaks=100, plot=F))
  new <- lapply(1:sample_num, function(x) hist(ps_var[,x], breaks=100, plot=F))
  newer <- lapply(1:sample_num, function(x) hist(ps_filt[,x], breaks=100, plot=F))
  plot(old[[x]], ylim=c(0,1000), xlab=samples[x], ylab = "", main = "", 
       border = "gray", col = "gray")
  plot(new[[x]], add=T, border = "green4", col = "green4")
  plot(newer[[x]], add=T, border = "red3", col = "green4")})
mtext("Pseudo Filtered with Medium Var", side = 3, line = -1, outer = T) 
p_pseudofvarhist <- recordPlot()
par(mfrow=c(1,1)) # Return to default page orientation
par(mar = c(4,4,4,4), oma = c(1,1,1,1))

# GSEA Filt ----
fac_DEseq <- grouped
coldata <- data.frame(cbind(fac_DEseq))
row.names(coldata) <- samples # DESeq2 uses raw counts; rows:genes & cols:samples
dds <- DESeqDataSetFromMatrix(countData=r_filt, colData=coldata, design = ~ fac_DEseq)
paste(nrow(dds), " genes input into DESeq2 for GSEA", sep="")
dds <- DESeq(dds)
human_genes1 <- convertMouseGeneList(as.array(rownames(get_counts(dds))))
human_genes1 <- human_genes1[!duplicated(human_genes1[1]),]
human_genes2 <- as.data.frame(human_genes1[,2])
row.names(human_genes2) <- human_genes1[,1]
colnames(human_genes2) <- "human_genes"
combined_list <- merge(get_counts(dds), human_genes2, by = "row.names", all.x = T)
nulllist <- combined_list[is.na(combined_list$human_genes),]
nulllist$human_genes <- nulllist$Row.names
notnulllist <- combined_list[!is.na(combined_list$human_genes),]
combined_list2 <- rbind(nulllist, notnulllist)
gct_format <- cbind(combined_list2[,length(combined_list2)], NA, 
                    combined_list2[,3:length(combined_list2)-1])
gct_format <- rbind("", "", colnames(gct_format), gct_format)
gct_format[1,1] <- "#1.2"
gct_format[2,1] <- nrow(combined_list2)
gct_format[2,2] <- ncol(combined_list2)-2
gct_format[3,1] <- "NAME"
gct_format[3,2] <- "Description"
write.table(gct_format, file = paste(currentDate, sample_set, "_filt.gct", sep=""), 
            quote=F, row.names=F, col.names=F, sep ="\t")
cls_format <- rbind("", "", colnames(combined_list2[3:ncol(combined_list2)-1]))
cls_format[3,] <- fac_DEseq
groups1 <- fac_DEseq[!duplicated(fac_DEseq)]
cls_format[1,1:3] <- c(ncol(cls_format), length(groups1),1)
cls_format[2,1:(length(groups1)+1)] <- c("#", groups1)
write.table(cls_format, file = paste(currentDate,sample_set, "_filt_classes.cls", 
                                     sep=""), quote=F, row.names=F, col.names=F)

# C) DIFF EXPRESSION ----------------------------------------------------------
# 9 T tests, Filt ----
# Group Filt Samples
Ad00 <- ps_filt[ , (cohort == "AlloDay00")]
Ad07 <- ps_filt[ , (cohort == "AlloDay07")]
Ad13 <- ps_filt[ , (cohort == "AlloDay13")]
AAll <- ps_filt[ , (donor == "Allogeneic")]
Sd00 <- ps_filt[ , (cohort == "SynDay00")]
Sd07 <- ps_filt[ , (cohort == "SynDay07")]
Sd13 <- ps_filt[ , (cohort == "SynDay13")]
SAll <- ps_filt[ , (donor == "Syngeneic")]
Bothd00 <- ps_filt[ , (time == "Day00")]
Bothd07 <- ps_filt[ , (time == "Day07")]
Bothd13 <- ps_filt[ , (time == "Day13")]
Health <- cbind(SAll, Ad00)
Reject <- cbind(Ad07, Ad13)
Tcells <- ps_filt[ , (immune == "Immunized")]
pairs <- list(Ad00, Ad07, Ad13, Sd00, Sd07, Sd13, AAll, SAll, Bothd00, Bothd07, 
              Bothd13, Health, Reject, Tcells)
pairs_names <- c("Ad00", "Ad07", "Ad13", "Sd00", "Sd07", "Sd13", "AAll", "SAll", 
                 "Bothd00", "Bothd07", "Bothd13", "Health", "Reject", "Tcells")
# Deltas Groups
A7min0 <- Ad07 - Ad00
colnames(A7min0) <- c("A7-0m5", "A7-0m6", "A7-0m8")
A13min0 <- Ad13 - Ad00
colnames(A13min0) <- c("A13-0m5", "A13-0m6", "A13-0m8")
A13min7 <- Ad13 - Ad07
colnames(A13min7) <- c("A13-7m5", "A13-7m6", "A13-7m8")
S7min0 <- Sd07 - Sd00
colnames(S7min0) <- c("S7-0m5", "S7-0m6", "S7-0m8")
S13min0 <- Sd13 - Sd00
colnames(S13min0) <- c("S13-0m5", "S13-0m6", "S13-0m8")
S13min7 <- Sd13 - Sd07
colnames(S13min7) <- c("S13-7m5", "S13-7m6", "S13-7m8")
pairsd <- list(A7min0, A13min0, A13min7, S7min0, S13min0, S13min7)
pairsd_names <-c("A7min0", "A13min0", "A13min7", "S7min0", "S13min0", "S13min7")
delta_filt <- cbind(S7min0, S13min0, S13min7, A7min0, A13min0, A13min7)
samples_delta <- colnames(delta_filt)
grouped_delta <- c(rep("Healthy", 9), rep("Early_Rejection", 3), 
                   rep("Late_Rejection", 3), rep("Progressing Rejection", 3))
simple_delta <- c(rep("Healthy", 9), rep("Rejecting", 9))
# Very cool stuff! But some of these custom functions I've made are VERY SLOW!
  # P values, All Groups - VERY SLOW! (15 min each)
  sig_g <- sig_pvaltable(pairs, pairs_names)
  write.csv(sig_g, file=paste(currentDate, filt_low, ",", filt_high, "filt",
                              "_all_t-test_p-value.csv", sep=""))
  sig_gdelta <- sig_pvaltable(pairsd, pairsd_names)
  write.csv(sig_gdelta, file=paste(currentDate, filt_low, ",", filt_high, "filt",
                                   "_all_t-test_p-value_delta.csv", sep=""))
par(mfrow=c(3,2))
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
A0to7 <- pvalgenes(Ad00, Ad07, plot=T)
A0to13 <- pvalgenes(Ad00, Ad13, plot=T)
A7to13 <- pvalgenes(Ad07, Ad13, plot=T)
A0toR <- pvalgenes(Ad00, Reject, plot=T)
A0toS0 <- pvalgenes(Ad00, Sd00, plot=T)
A7toS7 <- pvalgenes(Ad07, Sd07, plot=T)
mtext("Significant Genes, filtered CPM", side = 1, line = -1, outer = T)
p_pval1 <- recordPlot()
A13toS13 <- pvalgenes(Ad13, Sd13, plot=T)
AtoS <- pvalgenes(AAll, SAll, plot=T)
HtoR <- pvalgenes(Health, Reject, plot=T)
HtoA7 <- pvalgenes(Health, Ad07, plot=T)
HtoA13 <- pvalgenes(Health, Ad13, plot=T)
A7toB0 <- pvalgenes(Ad07, Bothd00, plot=T)
mtext("Significant Genes, filtered CPM", side = 1, line = -1, outer = T)
p_pval2 <- recordPlot()
A7d0toA13d0 <- pvalgenes(A7min0, A13min0, plot=T)
A7d0toA13d7 <- pvalgenes(A7min0, A13min7, plot=T)
A13d7toA13d0 <- pvalgenes(A13min7, A13min0, plot=T)
A7d0toS7d0 <- pvalgenes(A7min0, S7min0, plot=T)
A13d0toS13d0 <- pvalgenes(A13min0, S13min0, plot=T)
A13d7toS13d7 <- pvalgenes(A13min7, S13min7, plot=T)
mtext("Significant Genes, filtered CPM", side = 1, line = -1, outer = T)
p_pval3 <- recordPlot()
par(mfrow=c(1,1))
S0to7 <- pvalgenes(Sd00, Sd07)
S0to13 <- pvalgenes(Sd00, Sd13)
S7to13 <- pvalgenes(Sd07, Sd13)
S7toB0 <- pvalgenes(Sd07, Bothd00)
S13toB0 <- pvalgenes(Sd13, Bothd00)
A0toS7 <- pvalgenes(Ad00, Sd07)
A0toS13 <- pvalgenes(Ad00, Sd13)
A0toS <- pvalgenes(Ad00, SAll)
A7toS0 <- pvalgenes(Ad07, Sd00)
A7toS13 <- pvalgenes(Ad07, Sd13)
A7toS <- pvalgenes(Ad07, SAll)
A13toS <- pvalgenes(Ad13, SAll)
A13toS0 <- pvalgenes(Ad13, Sd00)
A13toS7 <- pvalgenes(Ad13, Sd07)
A13toB0 <- pvalgenes(Ad13, Bothd00)
RtoS0 <- pvalgenes(Reject, Sd00)
RtoS7 <- pvalgenes(Reject, Sd07)
RtoS13 <- pvalgenes(Reject, Sd13)
RtoS <- pvalgenes(Reject, SAll)
RtoB0 <- pvalgenes(Reject, Bothd00)
TctoA0 <- pvalgenes(Tcells, Ad00)
TctoS0 <- pvalgenes(Tcells, Sd00)
TctoB0 <- pvalgenes(Tcells, Bothd00)
# PCA Plots, from interesting groups above
par(mfrow=c(2,2))
plot_pca(n_filt[A0to7,])
plot_pca(n_filt[A0to13,])
plot_pca(n_filt[A7to13,])
plot_pca(n_filt[A0toR,])
mtext("Filtered Norm", side = 1, line = -1, outer = T)
p_pca1 <- recordPlot()    
plot_pca(n_filt[AtoS,])
plot_pca(n_filt[HtoR,])
plot_pca(n_filt[HtoA7,])
plot_pca(n_filt[HtoA13,])
mtext("Filtered Norm", side = 1, line = -1, outer = T)
p_pca2 <- recordPlot()
plot_pca(n_filt[A7d0toA13d0,])
plot_pca(n_filt[A7d0toA13d7,])
plot_pca(n_filt[A13d7toA13d0,])
plot_pca(n_filt[A7d0toS7d0,])
mtext("Filtered Norm", side = 1, line = -1, outer = T)
p_pca3 <- recordPlot()
plot_pca(n_filt[A13d0toS13d0,])
plot_pca(n_filt[A13d7toS13d7,])
mtext("Filtered Norm", side = 1, line = -1, outer = T)
p_pca4 <- recordPlot()  
# Venn Diagrams
par(mfrow=c(2,3))
venn(list("A0to7"=A0to7, "A0to13"=A0to13, "A7to13"=A7to13))
title("Allo to Allo Genes")
i1 <- Reduce(intersect, list(A0to7, A0to13, A7to13))
venn(list("A7toS0"=A7toS0, "A7toS7"=A7toS7, "A7toS13"=A7toS13))
title("Allo Day 7 to Syn Genes")
i2 <- Reduce(intersect, list(A7toS0, A7toS7, A7toS13))
venn(list("A13toS0"=A13toS0, "A13toS7"=A13toS7, "A13toS13"=A13toS13))
title("Allo Day 13 to Syn Genes")
i3 <- Reduce(intersect, list(A13toS0, A13toS7, A13toS13))
plot_pca(n_filt[Reduce(intersect, list(A0to7, A0to13, A7to13)), ])
plot_pca(n_filt[Reduce(intersect, list(A7toS0, A7toS7, A7toS13)), ])
plot_pca(n_filt[Reduce(intersect, list(A13toS0, A13toS7, A13toS13)), ])
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
heatmap.2(t(t(n_filt[i1, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_filt[i2, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_filt[i3, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
venn(list("HtoA7"=HtoA7, "HtoA13"=HtoA13, "HtoR"=HtoR))
title("All Healthy to Allo Genes")
i4 <- Reduce(intersect, list(HtoA7, HtoA13, HtoR))
venn(list("A7toB0"=A7toB0, "A13toB0"=A13toB0, "RtoB0"=RtoB0))
title("Rejecting to Day 0 Genes")
i5 <- Reduce(intersect, list(A7toB0, A13toB0, RtoB0))
venn(list("RtoS0"=RtoS0, "RtoS7"=RtoS7, "RtoS13"=RtoS13, "RtoS"=RtoS))
title("Rejecting to Syn Genes")
i6 <- Reduce(intersect, list(RtoS0, RtoS7, RtoS13, RtoS))
plot_pca(n_filt[Reduce(intersect, list(HtoA7, HtoA13, HtoR)), ])
plot_pca(n_filt[Reduce(intersect, list(A7toB0, A13toB0, RtoB0)), ])
plot_pca(n_filt[Reduce(intersect, list(RtoS0, RtoS7, RtoS13, RtoS)), ])
heatmap.2(t(t(n_filt[i4, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_filt[i5, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_filt[i6, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
par(mfrow=c(2,2))
venn(list("S0to7"=S0to7, "S0to13"=S0to13, "S7to13"=S7to13))
title("Syn to Syn Genes")
e1 <- Reduce(intersect, list(S0to7, S0to13, S7to13))
venn(list("A0toS0"=A0toS0, "A0toS7"=A0toS7, "A0toS13"=A0toS13, 
          "A0toS"=A0toS))
title("Allo Day 0 to Syn Genes")
e2 <- Reduce(intersect, list(A0toS0, A0toS7, A0toS13, A0toS))
plot_pca(n_filt[Reduce(intersect, list(S0to7, S0to13, S7to13)), ])
plot_pca(n_filt[Reduce(intersect, list(A0toS0, A0toS7, A0toS13, A0toS)), ])
heatmap.2(t(t(n_filt[e1, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_filt[e2, ])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
# Additional Filtering
include <- unique(c(A0to7, A0to13, A7to13, A7toS0, A7toS7, A7toS13, A13toS0, A13toS7, 
            A13toS13, HtoA7, HtoA13, A7toB0, A13toB0, RtoS0, RtoS7, RtoS13, RtoS, 
            RtoB0, HtoR, A7toS, A13toS))
exclude <- unique(c(S0to7, S0to13, S7to13, A0toS0, A0toS7, A0toS13, A0toS, S7toB0, 
                    S13toB0, TctoA0, TctoS0, TctoB0))
length(setdiff(include, exclude))
n_ttest <- n_counts[setdiff(include, exclude),]
g_ttest <- rownames(n_ttest) # Vector of filtered genes
write.csv(g_ttest, file= paste(currentDate, "_ttest_STx.csv"))
r_ttest <- r_counts[g_ttest,] # Used for DESeq2, VERY IMPORTANT
ps_ttest <- ps_counts[g_ttest,]
delta_ttest <- delta_filt[g_ttest,]
par(mfrow=c(1,1))
p_addfPCA <- plot.pca(n_ttest, grouped, "CPM interesting genes", ellipses=T)
ggbiplot(prcomp(t(n_ttest), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("Pseudo T-test") + theme_classic()

# GSEA T tests ----
# fac_DEseq <- grouped
# coldata <- data.frame(cbind(fac_DEseq))
# row.names(coldata) <- samples # DESeq2 uses raw counts; rows:genes & cols:samples
# dds <- DESeqDataSetFromMatrix(countData=r_ttest, colData=coldata, design = ~ fac_DEseq)
# paste(nrow(dds), " genes input into DESeq2 for GSEA", sep="")
# dds <- DESeq(dds)
# human_genes1 <- convertMouseGeneList(as.array(rownames(get_counts(dds))))
# human_genes1 <- human_genes1[!duplicated(human_genes1[1]),]
# human_genes2 <- as.data.frame(human_genes1[,2])
# row.names(human_genes2) <- human_genes1[,1]
# colnames(human_genes2) <- "human_genes"
# combined_list <- merge(get_counts(dds), human_genes2, by = "row.names", all.x = T)
# nulllist <- combined_list[is.na(combined_list$human_genes),]
# nulllist$human_genes <- nulllist$Row.names
# notnulllist <- combined_list[!is.na(combined_list$human_genes),]
# combined_list2 <- rbind(nulllist, notnulllist)
# gct_format <- cbind(combined_list2[,length(combined_list2)], NA, 
#                     combined_list2[,3:length(combined_list2)-1])
# gct_format <- rbind("", "", colnames(gct_format), gct_format)
# gct_format[1,1] <- "#1.2"
# gct_format[2,1] <- nrow(combined_list2)
# gct_format[2,2] <- ncol(combined_list2)-2
# gct_format[3,1] <- "NAME"
# gct_format[3,2] <- "Description"
# write.table(gct_format, file = paste(currentDate, sample_set, "_ttest.gct", sep=""), 
#             quote=F, row.names=F, col.names=F, sep ="\t")
# cls_format <- rbind("", "", colnames(combined_list2[3:ncol(combined_list2)-1]))
# cls_format[3,] <- fac_DEseq
# groups1 <- fac_DEseq[!duplicated(fac_DEseq)]
# cls_format[1,1:3] <- c(ncol(cls_format), length(groups1),1)
# cls_format[2,1:(length(groups1)+1)] <- c("#", groups1)
# write.table(cls_format, file = paste(currentDate,sample_set, "_ttest_classes.cls", 
#                                      sep=""), quote=F, row.names=F, col.names=F)

# 10 GSEA Results ----
# for the following code to work, you need to view the GSEA html index and then save the 
# tsv files. The negative tsv file rows are then manually copied over to the positive file.
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_16_Filt_All"
# name_dataset <- "STx_Filt_All"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Filt_Hallmark"
# name_dataset <- "STx_Filt_Hallmark"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Filt_Reactome"
# name_dataset <- "STx_Filt_Reactome"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Filt_GO"
# name_dataset <- "STx_Filt_GO"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Filt_ImmSign"
# name_dataset <- "STx_Filt_ImmSign"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Filt_CellType"
# name_dataset <- "STx_Filt_CellType"
folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_16_Ttest_All"
name_dataset <- "STx_Ttest_All"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Ttest_Hallmark"
# name_dataset <- "STx_Ttest_Hallmark"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Ttest_Reactome"
# name_dataset <- "STx_Ttest_Reactome"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Ttest_GO"
# name_dataset <- "STx_Ttest_GO"
# folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_08_Ttest_ImmSign"
# name_dataset <- "STx_Ttest_ImmSign"
num_genesets <- 400
# Load NES Scores
namefile1 <- paste(folder, "\\RvH.tsv", sep = "")
namefile2 <- paste(folder, "\\EvH.tsv", sep = "")
namefile3 <- paste(folder, "\\LvH.tsv", sep = "")
namefile4 <- paste(folder, "\\EvL.tsv", sep = "")
# DFs of rows=gene sets
nes_table1 <- read.table(namefile1, sep = "\t", header = T, fill = T)
nes_table1 <- nes_table1[, -c(2,3,4,5,9,10,11,12)] # Remove unneeded info
nes_table1 <- nes_table1[order(-nes_table1$NES),] # Order by NES values

nes1 <- rbind(head(nes_table1, n = num_genesets),
              tail(nes_table1, n = num_genesets)) # Top NES gene sets
nes_table2 <- read.table(namefile2, sep = "\t", header = T, fill = T)
nes_table2 <- nes_table2[, -c(2,3,4,5,9,10,11,12)]
nes_table2 <- nes_table2[order(-nes_table2$NES),]

nes2 <- rbind(head(nes_table2, n = num_genesets), 
              tail(nes_table2, n = num_genesets))
nes_table3 <- read.table(namefile3, sep = "\t", header = T, fill = T)
nes_table3 <- nes_table3[, -c(2,3,4,5,9,10,11,12)]
nes_table3 <- nes_table3[order(-nes_table3$NES),]

nes3 <- rbind(head(nes_table3, n = num_genesets), 
              tail(nes_table3, n = num_genesets))
nes_table4 <- read.table(namefile4, sep = "\t", header = T, fill = T)
nes_table4 <- nes_table4[, -c(2,3,4,5,9,10,11,12)]
nes_table4 <- nes_table4[order(-nes_table4$NES),]

nes4 <- rbind(head(nes_table4, n = num_genesets), 
              tail(nes_table4, n = num_genesets))
gs1 <- nes_table1[,1]
gs2 <- nes_table2[,1]
gs3 <- nes_table3[,1]
gs4 <- nes_table4[,1]
genesets1 <- nes1[,1]
genesets2 <- nes2[,1]
genesets3 <- nes3[,1]
genesets4 <- nes4[,1]
nes1[,1] <- genesets1
nes2[,1] <- genesets2
nes3[,1] <- genesets3
nes4[,1] <- genesets4
rownames(nes_table1) <- gs1
rownames(nes_table2) <- gs2
rownames(nes_table3) <- gs3
rownames(nes_table4) <- gs4
genesets <- unique(c(genesets1, genesets2, genesets3, genesets4))
tab1 <- nes_table1[genesets,]
tab2 <- nes_table2[genesets,]
tab3 <- nes_table3[genesets,]
tab4 <- nes_table4[genesets,]
nes_top <- cbind.data.frame(tab1, tab2[,2:4], tab3[,2:4], tab4[,2:4])
colnames(nes_top) <- c("GENESET", "Rejecting v Healthy", "p1", "q1", 
                       "Early Rejecting v Healthy", "p2", "q2", 
                       "Late Rejecting v Healthy", "p3", "q3", 
                       "Early Rejecting v Late", "p4", "q4")
write.csv(nes_top, file = paste(name_dataset, ".csv", sep = ""), row.names = F)
genesets <- substr(genesets, 1, 60) # Extract first 60 characters from gene set name
nes_top[,1] <- genesets
rownames(nes_top) <- genesets
nes_p <- nes_top[,-c(2,4,5,7,8,10,11,13)]
nes_pLF <- melt(nes_p, id = c("GENESET"))
nes_top <- nes_top[,-c(3,4,6,7,9,10,12,13)]
nes_topLF <- melt(nes_top, id = c("GENESET"))
nes_topLF$variable <- factor(nes_topLF$variable,levels=unique(nes_topLF$variable))
nes_topLF <- cbind.data.frame(nes_topLF, nes_pLF[,3])
colnames(nes_topLF) <- c("GENESET", "variable", "NES", "pval")
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57")
xx <- ggplot(nes_topLF, aes(x = GENESET, y = variable, label=round(NES, digits=2))) + 
  geom_point(aes(size = pval, fill = NES), alpha = 0.9, shape = 21) + 
  scale_fill_fermenter(palette = "RdBu", n.breaks = 10) +
  scale_size(trans = "reverse") +
  geom_text(size=3, nudge_y = 0.35, 
            color = ifelse(nes_topLF$NES > 1.50, 2, ifelse(nes_topLF$NES < -1.50, 4, 1))) + 
  scale_color_manual(values = c("deepred", "black")) + 
  labs( x= "Geneset", y = "Comparison", size = "P Value", fill = "NES")  + 
  # theme_light() +
  theme(axis.text.x = element_text(angle=30, hjust=0.95,vjust=1.05)) + coord_flip() + 
  theme(panel.grid.major.y = element_blank(), panel.border = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
xx
# After doing this, go into the csv file that was created, convert it to an excel file, and
# identify the highest/most interesting NES scores and genesets. I searched for note-
# worthy sets with NES >+-1.5 using these search terms: skin, allo, graft, transplant, antigen,
# immune, adapt, innate, cytokine, mono, myeloid, macrophage, neutrophil, nk, dc, dendritic,
# tcell, treg, tconv, spleen, apc, mdsc, reactome, kegg, gobp, hallmark, kine, thymus, LN, ly,
# lymph
# Next perform a Leading Edge Analysis, save the gmx file and convert it to txt:
lead <- 
  read.table(".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\leading_edge_matrix_154Sets_Filt_All_RvH.txt", 
             sep ="\t", header=T) # Leading edge .gmx file
lead <- lead[-1,] # Remove na row for each geneset
table(unlist(lead))
occurences <- table(unlist(lead))
occurences <- occurences[-1]
occurences <- occurences[order(-as.numeric(occurences))]
names(occurences)
sum(occurences > 2)
occurences[1:sum(occurences > 2)]
genes_lead <- names(occurences[1:sum(occurences > 2)])
data_frame <- filter(combined_list, human_genes %in% genes_lead) 
genes_leadm <- data_frame[,1] # Match human leading edge genes to their mouse counterparts
venn(list("Leading Edge" = genes_leadm, "T tests" = g_ttest))
title("Leading Edge Overlap")
genes_leadt <- unique(c(genes_leadm, g_ttest))

# Put Blue Genesets back in for nice image
folder <- ".\\Inputs\\STX_TC\\Mar15thGSEAAnalysis\\03_16_Ttest_All"
name_dataset <- "STx_Ttest_All"

num_genesets <- 20
# Load NES Scores
namefile1 <- paste(folder, "\\RvH.tsv", sep = "")
# namefile2 <- paste(folder, "\\EvH.tsv", sep = "")
# namefile3 <- paste(folder, "\\LvH.tsv", sep = "")
# namefile4 <- paste(folder, "\\EvL.tsv", sep = "")
# DFs of rows=gene sets
nes_table1 <- read.table(namefile1, sep = "\t", header = T, fill = T)
nes_table1 <- nes_table1[, -c(2,3,4,5,9,10,11,12)] # Remove unneeded info
nes_table1 <- nes_table1[order(-nes_table1$NES),] # Order by NES values

nes1 <- rbind(head(nes_table1, n = num_genesets),
              tail(nes_table1, n = num_genesets)) # Top NES gene sets
nes_table2 <- read.table(namefile2, sep = "\t", header = T, fill = T)
nes_table2 <- nes_table2[, -c(2,3,4,5,9,10,11,12)]
nes_table2 <- nes_table2[order(-nes_table2$NES),]

nes2 <- rbind(head(nes_table2, n = num_genesets), 
              tail(nes_table2, n = num_genesets))
nes_table3 <- read.table(namefile3, sep = "\t", header = T, fill = T)
nes_table3 <- nes_table3[, -c(2,3,4,5,9,10,11,12)]
nes_table3 <- nes_table3[order(-nes_table3$NES),]

nes3 <- rbind(head(nes_table3, n = num_genesets), 
              tail(nes_table3, n = num_genesets))
nes_table4 <- read.table(namefile4, sep = "\t", header = T, fill = T)
nes_table4 <- nes_table4[, -c(2,3,4,5,9,10,11,12)]
nes_table4 <- nes_table4[order(-nes_table4$NES),]

nes4 <- rbind(head(nes_table4, n = num_genesets), 
              tail(nes_table4, n = num_genesets))
gs1 <- nes_table1[,1]
gs2 <- nes_table2[,1]
gs3 <- nes_table3[,1]
gs4 <- nes_table4[,1]
genesets1 <- nes1[,1]
genesets2 <- nes2[,1]
genesets3 <- nes3[,1]
genesets4 <- nes4[,1]
nes1[,1] <- genesets1
nes2[,1] <- genesets2
nes3[,1] <- genesets3
nes4[,1] <- genesets4
rownames(nes_table1) <- gs1
rownames(nes_table2) <- gs2
rownames(nes_table3) <- gs3
rownames(nes_table4) <- gs4
genesets <- unique(c(genesets1, genesets2, genesets3, genesets4))
tab1 <- nes_table1[genesets,]
tab2 <- nes_table2[genesets,]
tab3 <- nes_table3[genesets,]
tab4 <- nes_table4[genesets,]
nes_top <- cbind.data.frame(tab1, tab2[,2:4], tab3[,2:4], tab4[,2:4])
colnames(nes_top) <- c("GENESET", "Rejecting v Healthy", "p1", "q1", 
                       "Early Rejecting v Healthy", "p2", "q2", 
                       "Late Rejecting v Healthy", "p3", "q3", 
                       "Early Rejecting v Late", "p4", "q4")
write.csv(nes_top, file = paste(name_dataset, ".csv", sep = ""), row.names = F)
genesets <- substr(genesets, 1, 60) # Extract first 60 characters from gene set name
nes_top[,1] <- genesets
rownames(nes_top) <- genesets
nes_p <- nes_top[,-c(2,4,5,7,8,10,11,13)]
nes_pLF <- melt(nes_p, id = c("GENESET"))
nes_top <- nes_top[,-c(3,4,6,7,9,10,12,13)]
nes_topLF <- melt(nes_top, id = c("GENESET"))
nes_topLF$variable <- factor(nes_topLF$variable,levels=unique(nes_topLF$variable))
nes_topLF <- cbind.data.frame(nes_topLF, nes_pLF[,3])
colnames(nes_topLF) <- c("GENESET", "variable", "NES", "pval")
colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57")
xx <- ggplot(nes_topLF, aes(x = GENESET, y = variable, label=round(NES, digits=2))) + 
  geom_point(aes(size = pval, fill = NES), alpha = 0.9, shape = 21) + 
  scale_fill_fermenter(palette = "RdBu", n.breaks = 10) +
  scale_size(trans = "reverse") +
  geom_text(size=3, nudge_y = 0.35, 
            color = ifelse(nes_topLF$NES > 1.50, 2, ifelse(nes_topLF$NES < -1.50, 4, 1))) + 
  scale_color_manual(values = c("deepred", "black")) + 
  labs( x= "Geneset", y = "Comparison", size = "P Value", fill = "NES")  + 
  # theme_light() +
  theme(axis.text.x = element_text(angle=30, hjust=0.95,vjust=1.05)) + coord_flip() + 
  theme(panel.grid.major.y = element_blank(), panel.border = element_blank(),
        # axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
xx

# 11 Elastic Net ----
# Sources
    # https://jasminedaly.com/tech-short-papers/glmnet_lasso_tutorial.html
    # https://www-frontiersin-org.proxy.lib.umich.edu/articles/10.3389/fgene.2013.00270/full#h3
    # https://github.com/KlinkeLab/ImmClass2019/blob/master/3.Classifier.R
# Define Factors
groupedf <-factor(grouped, labels = c(1,0,2)) # Factor(s) for Elastic Net (EN) later
simplef <-factor(simple, labels = c(0,1))
xfactors <- cbind.data.frame(groupedf)
rownames(xfactors) <- samples
# Load the data
dataset <- t(r_counts[genes_leadt,]) # Should this be the raw or some kind of normalized data!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
# View(dataset)
# Split the data into training and test set
set.seed(123)
training.samples <- dataset$groupedf %>% createDataPartition(p = 1.0, list = F)
train.data  <- dataset[training.samples, ]
test.data <- dataset[-training.samples, ]
x <- model.matrix(groupedf~., train.data)[,-1] # Predictor variables
y <- train.data$groupedf # Outcome variable
# Finding lambdas and alphas https://asmquantmacro.com/2016/04/26/fitting-elastic-net-model-in-r/
lambdagrid <- seq(0, 100)
alphagrid <- seq(0,1, length = 10)
srchgrid <- expand.grid(.alpha = alphagrid, .lambda = lambdagrid)
###### srch <- train()
genes_EN <- multinom3_EN( x, y, 100)
#LOOCV is what we already did?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!??!?!?!?!?!?!????!?!!!!!!!!?!!?
g_EN1 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`1`])
g_EN.95 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.95`])
g_EN.9 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.9`])
g_EN.85 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.85`])
g_EN.8 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.8`])
g_EN.75 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.75`])
g_EN.7 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.7`])
g_EN.65 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.65`])
g_EN.6 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.6`])
g_EN.55 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.55`])
g_EN.5 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.5`])
g_EN.45 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.45`])
g_EN.4 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.4`])
g_EN.35 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.35`])
g_EN.3 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.3`])
g_EN.25 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.25`])
g_EN.2 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.2`])
g_EN.15 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.15`])
g_EN.1 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.1`])
g_EN.05 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0.0499999999999999`])
g_EN0 <- unique(rownames(r_ttest)[genes_EN$ENgenes$`0`])
# Heatmaps and PCA plots
n_EN1 <- na.omit(n_counts[g_EN1,])
n_EN.95 <- na.omit(n_counts[g_EN.95,])
n_EN.9 <- na.omit(n_counts[g_EN.9,])
n_EN.85 <- na.omit(n_counts[g_EN.85,])
n_EN.8 <- na.omit(n_counts[g_EN.8,])
n_EN.75 <- na.omit(n_counts[g_EN.75,])
n_EN.7 <- na.omit(n_counts[g_EN.7,])
n_EN.65 <- na.omit(n_counts[g_EN.65,])
n_EN.6 <- na.omit(n_counts[g_EN.6,])
n_EN.55 <- na.omit(n_counts[g_EN.55,])
n_EN.5 <- na.omit(n_counts[g_EN.5,])
n_EN.45 <- na.omit(n_counts[g_EN.45,])
n_EN.4 <- na.omit(n_counts[g_EN.4,])
n_EN.35 <- na.omit(n_counts[g_EN.35,])
n_EN.3 <- na.omit(n_counts[g_EN.3,])
n_EN.25 <- na.omit(n_counts[g_EN.25,])
n_EN.2 <- na.omit(n_counts[g_EN.2,])
n_EN.15 <- na.omit(n_counts[g_EN.15,])
n_EN.1 <- na.omit(n_counts[g_EN.1,])
n_EN.05 <- na.omit(n_counts[g_EN.05,])
n_EN0 <- na.omit(n_counts[g_EN0,])
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(25)
colSide <- grouped
for (x in 1:sample_num) {
  if (grouped[x] == "Healthy") {
    colSide[x] <- "blue"
  } else if (grouped[x] == "Early_Rejecting") {
    colSide[x] <- "deeppink2"
  } else {
    colSide[x] <- "darkred" # Dark Red is Late rejecting
  }  }
heatmap.2(as.matrix(n_EN1), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx LASSO Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN1), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx LASSO Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.9), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.9 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.9), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.9 Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.8), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.8 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.8), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.8 Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.7), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.7 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.7), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.7 Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.65), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.65 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.65), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.65 Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.6), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.6 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.6), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.6 Geneset") + theme_classic()
heatmap.2(as.matrix(n_EN.5), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx EN alpha=0.5 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(n_EN.5), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx EN alpha=0.5 Geneset") + theme_classic()
length(g_EN.8)
write.csv(g_EN.8, file=paste(currentDate,"_ElasticNetGenes.csv",sep=""))

# 12 Delta Elastic Net ----
# Sources
# https://jasminedaly.com/tech-short-papers/glmnet_lasso_tutorial.html
# https://www-frontiersin-org.proxy.lib.umich.edu/articles/10.3389/fgene.2013.00270/full#h3
# https://github.com/KlinkeLab/ImmClass2019/blob/master/3.Classifier.R
# Define Factors
groupeddelf <-factor(grouped_delta, labels = c(0,1,2,3)) # Factor(s) for Elastic Net (EN) later
simpledelf <-factor(simple_delta, labels = c(0,1))
xfactors <- cbind.data.frame(simpledelf)
rownames(xfactors) <- samples_delta
# Load the data
dataset <- t(delta_ttest) # Should this be the raw or some kind of normalized data!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
dataset <- as.data.frame(cbind(xfactors, dataset)) # Add factor as the first column
dataset <- dataset[complete.cases(dataset), ] # For cases without missed items
# View(dataset)
# Split the data into training and test set
set.seed(123)
training.samples <- dataset$simpledelf %>% createDataPartition(p = 1.0, list = F)
train.data  <- dataset[training.samples, ]
test.data <- dataset[-training.samples, ]
x <- model.matrix(simpledelf~., train.data)[,-1] # Predictor variables
y <- train.data$simpledelf # Outcome variable
# Finding lambdas and alphas https://asmquantmacro.com/2016/04/26/fitting-elastic-net-model-in-r/
lambdagrid <- seq(0, 100)
alphagrid <- seq(0,1, length = 10)
srchgrid <- expand.grid(.alpha = alphagrid, .lambda = lambdagrid)
###### srch <- train()
genes_EN <- binom_EN( x, y, 100)
#LOOCV is what we already did?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!??!?!?!?!?!?!????!
g_EN1 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`1`])
g_EN.95 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.95`])
g_EN.9 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.9`])
g_EN.85 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.85`])
g_EN.8 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.8`])
g_EN.75 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.75`])
g_EN.7 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.7`])
g_EN.65 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.65`])
g_EN.6 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.6`])
g_EN.55 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.55`])
g_EN.5 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.5`])
g_EN.45 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.45`])
g_EN.4 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.4`])
g_EN.35 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.35`])
g_EN.3 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.3`])
g_EN.25 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.25`])
g_EN.2 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.2`])
g_EN.15 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.15`])
g_EN.1 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.1`])
g_EN.05 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0.0499999999999999`])
g_EN0 <- unique(rownames(delta_ttest)[genes_EN$ENgenes$`0`])
# Heatmaps and PCA plots
d_EN1 <- na.omit(delta_ttest[g_EN1,])
d_EN.95 <- na.omit(delta_ttest[g_EN.95,])
d_EN.9 <- na.omit(delta_ttest[g_EN.9,])
d_EN.85 <- na.omit(delta_ttest[g_EN.85,])
d_EN.8 <- na.omit(delta_ttest[g_EN.8,])
d_EN.75 <- na.omit(delta_ttest[g_EN.75,])
d_EN.7 <- na.omit(delta_ttest[g_EN.7,])
d_EN.65 <- na.omit(delta_ttest[g_EN.65,])
d_EN.6 <- na.omit(delta_ttest[g_EN.6,])
d_EN.55 <- na.omit(delta_ttest[g_EN.55,])
d_EN.5 <- na.omit(delta_ttest[g_EN.5,])
d_EN.45 <- na.omit(delta_ttest[g_EN.45,])
d_EN.4 <- na.omit(delta_ttest[g_EN.4,])
d_EN.35 <- na.omit(delta_ttest[g_EN.35,])
d_EN.3 <- na.omit(delta_ttest[g_EN.3,])
d_EN.25 <- na.omit(delta_ttest[g_EN.25,])
d_EN.2 <- na.omit(delta_ttest[g_EN.2,])
d_EN.15 <- na.omit(delta_ttest[g_EN.15,])
d_EN.1 <- na.omit(delta_ttest[g_EN.1,])
d_EN.05 <- na.omit(delta_ttest[g_EN.05,])
d_EN0 <- na.omit(delta_ttest[g_EN0,])
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(25)
colSide <- grouped_delta
for (x in 1:sample_num) {
  if (grouped_delta[x] == "Healthy") {
    colSide[x] <- "blue"
  } else if (grouped_delta[x] == "Early_Rejection") {
    colSide[x] <- "deeppink2"
  } else if (grouped_delta[x] == "Late_Rejection") {
    colSide[x] <- "red"
  } else {
    colSide[x] <- "darkred" # Dark Red is Late rejecting
  }  }
heatmap.2(as.matrix(d_EN1), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN LASSO Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN1), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN LASSO Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.9), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.9 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.9), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.9 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.8), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.8 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.8), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.8 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.7), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.7 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.7), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.7 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.65), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.65 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.65), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.65 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.6), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.6 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.6), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.6 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.5), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.5 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.5), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.5 Geneset") + theme_classic()

heatmap.2(as.matrix(d_EN.2), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.2 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.2), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.2 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.1), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.1 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.1), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.1 Geneset") + theme_classic()
heatmap.2(as.matrix(d_EN.05), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx delEN alpha=0.05 Geneset", # density.info="none",
          margins = c(7, 7), ColSideColors = colSide, trace="none", # Colv = "Na",
          key.title=NA, key.ylab=NA, keysize = 1.0, dendrogram = "column")
ggbiplot(prcomp(t(d_EN.05), scale.=T), ellipse=T, groups=grouped_delta, var.axes=F, 
         labels=samples_delta, var.scale=1, circle=T) + 
  ggtitle("STx delEN alpha=0.05 Geneset") + theme_classic()
length(g_EN.6)
write.csv(g_EN.6, file=paste(currentDate,"_DeltaElasticNetGenes.csv",sep=""))
par(mar=c(4,4,4,4))
par(mfrow=c(1,1))
venn(list("delEN 0.05" = rownames(d_EN.05), "EN 0.05" = rownames(n_EN.05)))
title("Elastic Net Overlap")
# Volcano plot
# genes_all <- unique(c(rownames(d_EN.05), rownames(n_EN.05)))
t.test(r_counts[rownames(d_EN.05),1:12], r_counts[rownames(d_EN.05),13:18])
ttest_all <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue <- apply(r_counts[rownames(d_EN.05),], 1, ttest_all, grp1 = c(1:12), grp2 = c(13:18))
hist(rawpvalue)
r_counts[rownames(d_EN.05),] <- log2(r_counts[rownames(d_EN.05),]) # Transform data into log2 base
control <- apply(r_counts[rownames(d_EN.05), 1:12], 1, mean) # Find mean of genes in control
test <- apply(r_counts[rownames(d_EN.05), 13:18], 1, mean) # Find mean of genes in test group
class(control) # Confirming that we have a vector of numbers
class(test)
foldchange <- control - test # bc data is log2 transformed, take difference btwn means
hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")
results = cbind(foldchange, rawpvalue) # log2FC or log2 Ratio=log2(control/test)
results = as.data.frame(results)
results$probename <- rownames(results)
volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point() + 
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=1.0) +
  geom_hline(yintercept=2, linetype="dashed", color = "darkred", size=1.0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -0.5, linetype="dotted", color = "blue", size=1.0) +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "blue", size=1.0) +
  geom_vline(xintercept = -0.75, linetype="dotted", color = "blue", size=1.0) +
  geom_vline(xintercept = 0.75, linetype="dotted", color = "darkblue", size=1.0)
# Color dots next  

# 13 DESeq2 ----
fac_DEseq <- grouped
coldata <- data.frame(cbind(fac_DEseq))
row.names(coldata) <- samples # DESeq2 uses raw counts; rows:genes & cols:samples
dds <- DESeqDataSetFromMatrix(countData=r_ttest, colData=coldata, design = ~ fac_DEseq)
paste(nrow(dds), " genes input into DESeq2 from T tests", sep="")
dds <- DESeq(dds)
res <- results(dds) # Table of log2 fold changes, p-values, & p-adj values
res # Inspect results tables
mcols(res, use.names = T) # View metadata of results
summary(res)
resultsordered <- res[order(res$padj),]   # Order results by p-adj
sum(res$padj < 0.1, na.rm=T)   # How many genes have p-adj < 0.1
resSig <- subset(resultsordered, padj < 0.1) # Subset by p-adj < 0.1
head(resSig[order(resSig$log2FoldChange),]) # Sig genes w/ strongest down-regulation
tail(resSig[order(resSig$log2FoldChange),]) # Sig genes w/ strongest up-regulation
dds_genes <- rownames(resSig) # Differential expression
r_DEG <- r_counts[dds_genes,] 
n_DEG <- n_counts[dds_genes,]
par(mfrow=c(1,1)) # Visualize Results
par(mar = c(4,4,4,4), oma = c(1,1,1,1)) # Adjust margins to better use space
plotMA(res,main=paste("Log2FC of DESEQ2,", toString(nrow(resSig)), 
                      "Genes padj<0.1,", filt_low, "&", filt_high, "filters"))
p_lfc <- recordPlot()
plotDispEsts(dds)
p_hist <- hist(res$pvalue, breaks=50, plot=F) # Store histogram
coul <- ifelse(p_hist$breaks<=0.05, "light green", "grey") # Colors
plot(p_hist, col=coul, border=F, main="", xlab="p value", ylab="Freq")
res$pvalue[is.na(res$pvalue)] <- 0 # Replace NaN's with zeros
sum(res$pvalue>0 & res$pvalue<0.05)
p_hist <- hist(res$padj, breaks=50, plot=F) # Store histogram
coul <- ifelse(p_hist$breaks<=0.1, "light green", "grey") # Colors
plot(p_hist, col=coul, border=F, main="", xlab="p adj", ylab="Freq")
res$padj[is.na(res$padj)] <- 0 # Replace NaN's with zeros
sum(res$padj>0 & res$padj<0.1)
qs <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:6/6))
bins <- cut(res$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(res$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")
LFC_cut <- 0.8
padj_cut <- 0.1
par(mar=c(4,4,2,1), cex=1.0, cex.main=1.4, cex.axis=1.0, cex.lab=1.2)
res_df <- as.data.frame(res)
g_cutoff <- rownames(subset(res_df, padj<padj_cut & padj>0 & abs(log2FoldChange)>LFC_cut))
with(res_df, plot( -log10(padj), log2FoldChange, pch=20, cex=1.0,
                   main=paste("Volcano plot, ", length(g_cutoff), " genes ±LFC>", 
                              LFC_cut, " & ", "padj<", padj_cut, sep=""), 
     xlab = bquote(~-log[10]~Q~value)), ylab = bquote(~Log[2]~fold~change))
with(subset(res_df, padj<padj_cut & abs(log2FoldChange)<LFC_cut), # padj < 0.1
     points( -log10(padj), log2FoldChange, pch=20, col="blue", cex=1))
with(subset(res_df, padj<padj_cut & abs(log2FoldChange)>LFC_cut), # LFC > 0.5
     points( -log10(padj), log2FoldChange, pch=20, col="orange3", cex=1))
with(subset(res_df, -log10(padj)>1.5 & abs(log2FoldChange)>LFC_cut), # LFC > 0.5
     points( -log10(padj), log2FoldChange, pch=20, col="red3", cex=1))
with(subset(res_df, padj<padj_cut & abs(log2FoldChange)>LFC_cut), # Label points
     text( -log10(padj), log2FoldChange,
          labels=subset(rownames(res_df), res_df$padj<padj_cut &
                          abs(res_df$log2FoldChange)>LFC_cut), cex=0.6, pos=3))
abline(h=0, col="black", lty=3, lwd=0.8) # Add center line
abline(h=-LFC_cut, col="blue", lty=2, lwd=1.5) # Line for -LFC cutoff
abline(h=LFC_cut, col="blue", lty=2, lwd=1.5) # Line for +LFC cutoff
abline(v=-log10(padj_cut), col="orange3", lty=2, lwd=1.5) # Line for p adj cutoff
abline(v=1.5, col="red3", lty=2, lwd=1.5) # Line for p adj cutoff
n_LFC <- n_counts[g_cutoff, ]
# Heatmap & PCA
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(25)
colSide <- grouped
for (x in 1:sample_num) {
  if (grouped[x] == "Healthy") {
    colSide[x] <- "blue"
  } else if (grouped[x] == "Early_Rejecting") {
    colSide[x] <- "deeppink2"
  } else {
    colSide[x] <- "darkred" # Dark Red is Late rejecting
  }  }
heatmap.2(as.matrix(n_LFC), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx LFC Geneset", 
          margins = c(7, 7), ColSideColors = colSide, trace="none",
          key.title=NA, key.ylab=NA, keysize = 1.0, # density.info="none",
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_LFC), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx LFC Geneset") + theme_classic()
g_topVar <- head(order(rowVars(n_LFC), decreasing = T), 12) # Top DEGs by variance
mat  <- n_LFC[ g_topVar, ]
mat  <- mat - rowMeans(mat)
heatmap.2(as.matrix(mat), scale = "row", col = coul, key = T, 
          xlab = "Scaffold", ylab="STx LFC Top 12 Geneset", 
          margins = c(7, 7), ColSideColors = colSide, trace="none",
          key.title=NA, key.ylab=NA, keysize = 1.0, # density.info="none",
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(mat), scale.=T), ellipse=T, groups=grouped, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("STx LFC Top 12 Geneset") + theme_classic()
# Looking at common genes
par(mfrow=c(1,1))
par(mar=c(2,2,2,2))
venn(list("Delta EN" = rownames(d_EN.5), "EN" = rownames(n_EN.5), 
          "DESep DEG" = rownames(n_DEG)))
title("Narrowed Down Genes")
genes_all <- unique(c(rownames(d_EN.05), rownames(n_EN.05), rownames(n_DEG)))

# 14 SVD ----
n_centered <- n_DEG - rowMeans(n_DEG) # First, center data on genes
svd1 <- svd(t(n_centered)) # Apply SVD on transposed data
Z <- t(n_centered) %*% svd1$v # Z=XV
svd1$u # U are unscaled PCs
par(mfrow=c(1,1)) # Only one plot per page
plot(svd1$u[,1], svd1$u[,2], col = coul, main = "Genes v Samples (SVA)", pch=19,
     xlab = "V1", ylab = "V2")
text(svd1$u[,1], svd1$u[,2], samples, cex=0.7, pos=4)
p_svd <- recordPlot()
centr_healthy <- colMeans(rbind(svd1$u[clin_info$state == 0,]),dims=1)
euclid <- rbind(centr_healthy, svd1$u) # Euclidean Distances
euclid_dist <- rdist(euclid[,1:3], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:18]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_SVD <- cbind(euc, mouse, time, cohort, donor, simple, grouped)
p_svdbar <- ggplot(euc, aes(x = simple, y = euclid_dist, fill = cohort)) + 
  geom_boxplot() + xlab("") + ylab("") + ggtitle("Singular Value Decomposition")
p_svdbar

# 15 Random Forest ----
predictor_data <- t(n_LFC) # Transpose data & assign genes as predictors.
target <- clin_info[,"state"] # Set variable to predict target (reject status)
target[target==0] <- "Healthy"
target[target==1] <- "Rejecting"
# target[target==1] <- "RejectingD7"
# target[target==2] <- "RejectingD13"
target <- as.factor(target)
# Odd # for ntree because ties are broken randomly, & an odd number of trees makes 
# the model fully deterministic. Down-sampling to compensate for unequal classes
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=F)[1]]
sampsizes <- rep(min_size,num_classes)
rf_output <- randomForest(x=predictor_data, y=target, importance=T, ntree=10001, 
                          proximity=T, sampsize=sampsizes, na.action=na.omit)
rf_importances <- importance(rf_output, scale=F) # Importance for each variable
confusion <- rf_output$confusion # Calculates sensitivity, specificity, accuracy
sensitivity <- (confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity <- (confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error <- rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy <- 1-overall_error
class1_error <- paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error <- paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy <- 100-overall_error
# Prepare Output File
sens_out <- paste("sensitivity=",sensitivity, sep="")
spec_out <- paste("specificity=",specificity, sep="")
err_out <- paste("overall error rate=",overall_error,sep="")
acc_out <- paste("overall accuracy=",overall_accuracy,sep="")
misclass_1 <- paste(confusion[1,2], rownames(confusion)[1],"misclassified as", 
                 colnames(confusion)[2], sep=" ")
misclass_2 <- paste(confusion[2,1], rownames(confusion)[2],"misclassified as", 
                 colnames(confusion)[1], sep=" ")
confusion_out <- confusion[1:2,1:2]
confusion_out <- cbind(rownames(confusion_out), confusion_out)
# Representation of top 30 variables, categorized by importance.
p_predictors <- varImpPlot(rf_output, type=2, n.var=12, scale=F, 
                    main="Variable Importance (Gini) for DEG predictors")
# MDS Class Separation
target_labels <- as.vector(target)
target_labels[target_labels=="Healthy"] <- "H"
target_labels[target_labels=="Rejecting"] <- "R"
plot_MDS <- MDSplot(rf_output, target,k=2,xlab="",ylab="",pch=target_labels, 
                    palette=c("red", "blue"), main="MDS plot")
plot_MDS <- plot_MDS$points
plot_MDS <- as.data.frame(plot_MDS)
p_mds <- ggplot(plot_MDS, aes(x=`Dim 1`,y=`Dim 2`, color=target_labels)) + 
  geom_point(aes(shape=donor),size=3) + 
  geom_text(aes(label = time),nudge_x=0.03, nudge_y=-0.01)
centr_healthy <- colMeans(rbind(plot_MDS[clin_info$state == 0,]),dims=1)
euclid <- rbind(centr_healthy, plot_MDS) # Euclidean distances
euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
euclid_dist <- euclid_dist[1:18]
euc <- cbind.data.frame(samples, euclid_dist)
rownames(euc) <- c(samples)
euc_RF <- cbind(euc, mouse, time, cohort, donor, simple, grouped)
p_rf <- ggplot(euc, aes(x = simple, y = euclid_dist, fill = cohort)) + 
  geom_boxplot() + xlab("") + ylab("") + ggtitle("Random Forest")
p_rf
ggarrange(plotlist=list(p_svdbar, p_rf), common.legend=T, 
          legend="right", nrow=1, ncol=2)
# ROC Curve, Rejecting/Healthy vote fractions as predictive variable. ROC curve 
# made by stepping through different thresholds for calling rejecting vs healthy.
predictions <- as.vector(rf_output$votes[,2])
pred <- prediction(predictions,target)
perf_AUC <- performance(pred,"auc") # First calculate the AUC value
AUC <- perf_AUC@y.values[[1]]
perf_ROC <- performance(pred,"tpr","fpr") # Plot the actual ROC curve
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=4, scientific=F)))
p_ROC <- recordPlot()
# Vote Distributions
options(digits=2)
out <- histbackback(split(rf_output$votes[,"Rejecting"], target), 
                    probability=F, axes=T, xlim=c(-50,50),
                    main='Vote distributions for mice classified by RF', 
                    ylab="Fraction votes (Rejecting)")
barplot(-out$left, col="red" , horiz=T, space=0, add=T, axes=F)
barplot(out$right, col="blue", horiz=T, space=0, add=T, axes=F)
p_votes <- recordPlot()

# 16 Scoring System ----
plot(euc_SVD$euclid_dist, euc_RF$euclid_dist)
scores <- cbind.data.frame(euc_SVD$euclid_dist, euc_RF$euclid_dist)
colnames(scores) <- c("SVD", "RF")
rownames(scores) <- samples
ggplot(scores, aes(x=SVD,y=RF, color=grouped)) + 
  geom_point(aes(shape=donor),size=3) + 
  geom_text(aes(label = time),nudge_x=0.03, nudge_y=-0.01) + 
  xlab("Singular Value Decomposition (arb. units)") + 
  ylab("Random Forest prediction, Rejecting probability") + 
  ggtitle(paste( nrow(n_LFC), " DEGs, Filtered CPM: ", filt_low, " low, ", 
                 filt_group, " group, & ", filt_high, " high filters from ", 
                 "SIMPLE", # CHANGE THIS AS NEEDED!
                 " factor", sep="")) + 
  # Add ggforce's ellipses
  ggforce::geom_mark_ellipse(aes(fill = grouped,color = grouped)) + coord_equal()
  # Need to add confidence intervals!
p_scoring <- recordPlot()
# Create Plot
coul <- c(rep("deepskyblue", 9), rep("darkorange3", 3), rep("deepskyblue", 3), 
          rep("red3", 3))
plot (scores$SVD, scores$RF, col=coul, pch=19, cex=2,
      ylab="Random Forest prediction (prob.)", 
      xlab="Singular Value Decomposition (arb. units)") 
#Add Confidence Interval Ellipses
dataEllipse(scores$SVD[clin_info$state == 0], scores$RF[clin_info$state == 0], 
            levels=c(0.7), center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points=FALSE, col="deepskyblue3", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[clin_info$grouped == "Early_Rejecting"], 
            scores$RF[clin_info$grouped == "Early_Rejecting"], levels=c(0.7), 
            center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, robust=FALSE, 
            plot.points = FALSE, col="darkorange4", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[clin_info$grouped == "Late_Rejecting"], 
            scores$RF[clin_info$grouped == "Late_Rejecting"], levels=c(0.7), 
            center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, robust=FALSE, 
            plot.points = FALSE, col="red4", pch=1, lwd=1, lty=1)
dataEllipse(scores$SVD[clin_info$state == 1], scores$RF[clin_info$state == 1], 
            levels=c(0.7), center.pch=FALSE, draw=TRUE, add=TRUE, segments=51, 
            robust=FALSE, plot.points=FALSE, col="red4", pch=1, lwd=1, lty=1)
# Add a legend
legend("bottomleft", 
       legend = unique(grouped), col = unique(coul), pch = 19, bty = "n", 
       pt.cex = 1, cex = 1, text.col = "black", horiz = F , 
       inset = c(0.005, 0.84))
# CSV of Signature
write.csv(rownames(n_LFC), file= "signature_STx.csv")

# 17 Time Course Exp ----
coldata <- data.frame(cbind(donor, day)) # Sample information
row.names(coldata) <- samples
ddsTC <- DESeqDataSetFromMatrix(countData = r_ttest, colData = coldata, 
                                design = ~ donor + day + donor:day)
nrow(ddsTC)
# Variables start with ~, with + signs between them. Put control variable first.
ddsTC <- DESeq(ddsTC, test = "LRT", reduced = ~ donor + day)
# Likelihood Ratio Test: remove donor response-specific differences over time.
# Genes with small p value, which after time 0 show a donor response-specific 
# effect. Doesn't give small p values to genes that changed over time in the 
# same way in both donor responses.
resTC <- results(ddsTC) # Table of log2 fold changes, p values, & p-adj values
resTC
mcols(resTC, use.names = T)
summary(resTC)
resTC$symbol <- mcols(ddsTC)$symbol # What does this line do?
head(resTC[order(resTC$padj),], 10)
# Organize Results
resTC_ordered <- resTC[order(res$padj),] # Order results by p-adj
sum(resTC$padj < 0.1, na.rm=T) # Subset by p-adj < 0.1
resTC_Sign <- subset(resTC_ordered, padj < 0.1)
ddsTC_genes <- rownames(resTC_Sign)
r_TCDEG <- r_counts[ddsTC_genes,]
n_TCDEG <- n_counts[ddsTC_genes,]
# Plot Results, counts over time for genes with low padj, accounting for condition-
# dependent time profile & time 0 differences. Interaction terms are difference 
# between donor groups at a given time after accounting for the difference at time 0.
plot_lst <- c() # Empty list of plots
for (i in 1:length(ddsTC_genes)) {
  time_course <- plotCounts(ddsTC, gene=ddsTC_genes[i], 
                            intgroup = c("day","donor"), returnData = T)
  time_course$day <- as.numeric(as.character(time_course$day))
  plot_lst[[i]] <- ggplot(time_course, aes(x=day, y=count, color=donor, 
                                           shape=donor, group=donor)) + 
    geom_point(size = 2) + 
    stat_summary(fun=mean, geom="line") +
    ggtitle(ddsTC_genes[i]) +
    scale_y_log10()
}
par(mar = c(4,1,2,1), oma = c(0.5,1,0,0)) # Adjust margins to better use space
ggarrange(plotlist=plot_lst, common.legend=T, legend="bottom", nrow=3, ncol=3)
time_data <- lapply(1:length(ddsTC_genes), function(x) 
  plotCounts(ddsTC, gene=ddsTC_genes[i], 
             intgroup = c("day","donor"), returnData = T))
names(time_data) <- ddsTC_genes
# time_data <- unlist(time_data)
resultsNames(ddsTC)
res30 <- results(ddsTC, name="donorSyngeneic.day7", test="Wald")
res30[which.min(resTC$padj),]
betas <- coef(ddsTC)
colnames(betas)
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks = seq(from = -thr, to = thr, length = 101), cluster_col = F)
# clusters <- degPatterns(cluster_rlog, metadata = meta, time="say", col="donor")
# Heatmap
coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
colSide <- c(rep("darkblue", 3), rep("darkred", 3), rep("deepskyblue", 3), 
             rep("orange", 3), rep("blue", 3), rep("deeppink2", 3))
heatmap.2(t(t(n_TCDEG)), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
p_TCheatmap <- recordPlot()
plot_pca(n_TCDEG) # PCA
p_TCpcaDEG <- recordPlot()

# Panther GO ----
panther_ttest <- read.table(".\\Inputs\\STX_TC\\addfilt_4-11-21.csv", sep =",", 
                          header=F, fill=T)
colnames(panther_ttest) <- c("Category",	"PANTHER GO-Slim Biological Process",	
                            "Genes",	"Percent of Genes",	"Percent of Hits")
ggplot(panther_ttest, aes(x=`PANTHER GO-Slim Biological Process`, y=Genes,
                         fill=`PANTHER GO-Slim Biological Process`)) + 
  geom_bar(stat="identity", width=1, color="white") +
  # theme(axis.text.x  = element_text(angle=30, hjust=0.95,vjust=1.05)) + 
  coord_flip() + 
  theme(
    legend.position = "none"
    # panel.grid.major.y = element_blank(),
    # panel.border = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    # axis.text.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank()
    )

panther_sign <- read.table(".\\Inputs\\STX_TC\\Signature_4-11-21.csv", sep =",", 
                           header=F, fill=T)
colnames(panther_sign) <- c("Category",	"Biological Process",	
                            "Genes",	"Percent of Genes",	"Percent of Hits")
ggplot(panther_sign, aes(x="", y=Genes,
                         fill=`PANTHER GO-Slim Biological Process`)) + 
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() # remove background, grid, numeric labels

ggplot(panther_sign, aes(x=`Biological Process`, y=Genes,
                            fill=`Biological Process`)) + 
  geom_bar(stat="identity", width=0.9, color="white") +
  theme(axis.text.y  = element_text(size=15), 
        axis.title.y  = element_text(size=18)) +
  coord_flip() + 
  labs( x= "Biological Process", y = "Genes") + 
  theme(
    legend.position = "none"
    # panel.grid.major.y = element_blank(),
    # panel.border = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    # axis.text.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.x = element_blank()
  )

# Housekeeping Genes ----
hk <- c("Cdkn1a", "Pgk1", "Tfrc", "Gusb", "Hprt1", "Tbp", "Hmbs", "Rplp2", "Polr2a", 
        "Ywhaz", "Ipo8", "Ppia", "Gapdh")
n_hk <- n_counts[hk,]
n_hk <- n_hk[rownames(n_hk) == hk, ]
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(25)
colSide <- grouped
for (x in 1:sample_num) {
  if (grouped[x] == "Healthy") {
    colSide[x] <- "blue"
  } else if (grouped[x] == "Early_Rejecting") {
    colSide[x] <- "deeppink2"
  } else {
    colSide[x] <- "darkred" # Dark Red is Late rejecting
  }  }
heatmap.2(t(t(n_hk)), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
n_probes <- rbind(n_DEG, n_hk)
heatmap.2(t(t(n_probes)), scale="row", col = coul, key=T, xlab="Mouse", 
        ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
        key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
par(mfrow=c(1,1)) # Return to default page orientation
var_cut <- 0.005
df <- cbind.data.frame("x"=rowMeans(ps_counts), "y"=rowVars(ps_counts))
g_lowvar <- rownames(subset(df, y < var_cut))
with(df, plot(x, y, pch=20, cex=1,
              main=paste("All Var vs Mean, ", length(g_lowvar), 
                         " genes with Var > ", var_cut, sep=""), 
              xlab="Mean across 18 samples", ylab="Var across 18 samples"))
with(subset(df, y < var_cut), 
     points(x, y, pch=20, col="red3", cex=1))
abline(h=var_cut, col="red3", lty=2, lwd=1.5) # Line for Var cutoff
ggbiplot(prcomp(t(ps_counts[ g_lowvar,]), scale.=T), ellipse=T, groups=cohort, 
         var.axes=F, labels=mouse, var.scale=1, circle=T) + 
  ggtitle("Pseudo: Lowest Var") + theme_classic() 
ggscreeplot(prcomp(t(ps_counts[ g_lowvar,]), scale.=T))
housekeep <- g_lowvar
housekeep <- Reduce(intersect, list(hk, g_lowvar))
housekeep <- c("Sypl", "Tnrc6a")
venn(list("OpenArray HK"=hk, "Low Var"=g_lowvar))
title("Finding HK Genes")
ggbiplot(prcomp(t(ps_counts[ housekeep,]), scale.=T), ellipse=T, groups=cohort, 
         var.axes=F, labels=mouse, var.scale=1, circle=T) + 
  ggtitle("Pseudo: Low Var") + theme_classic() 
ggscreeplot(prcomp(t(ps_counts[ g_lowvar,]), scale.=T))
heatmap.2(t(t(n_counts[housekeep,])), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
n_probes <- rbind(n_DEG, n_counts[housekeep,])
heatmap.2(t(t(n_probes)), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
heatmap.2(t(t(n_DEG)), scale="row", col = coul, key=T, xlab="Mouse", 
          ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
out <- heatmap.2(t(t(n_DEG)))
plot(out$rowDendrogram)
g_DEGorder <- rownames(mat)[out$rowInd]
s_DEGorder <- colnames(mat)[out$colInd]
g_probes <- c(g_DEGorder, housekeep)
heatmap.2(t(t(n_counts[g_probes, s_DEGorder])), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          Rowv = F, Colv = F, key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
mat <- n_counts[g_probes, s_DEGorder]
mat  <- mat - rowMeans(mat)
heatmap.2(t(t(mat)), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          Rowv = F, Colv = F, key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
oarray <- c("Hprt1", "Tbp", "Hmbs", "Rplp2", "Polr2a", "Ipo8", "Ppia", "Gapdh")
heatmap.2(t(t(n_counts[oarray, ])), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
g_probes <- c(g_DEGorder, oarray)
heatmap.2(t(t(n_counts[g_probes, s_DEGorder])), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          Rowv = F, Colv = F, key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
mat <- n_counts[g_probes, s_DEGorder]
mat  <- mat - rowMeans(mat)
heatmap.2(t(t(mat)), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          Rowv = F, Colv = F, key.title=NA, key.ylab=NA, trace="none", tracecol=NA)
hk <- c("Cdkn1a", "Pgk1", "B2m", "Tfrc", "Gusb", "Hprt1", "Tbp", "Actb", "Hmbs", 
        "Rplp2", "Polr2a", "Ywhaz", "Ubc", "Ipo8", "Ppia", "Gapdh")
ggplot(data=n_counts["Cdkn1a", ])
# Convert to "long form" table via melt function for generating some plots
colnames(n_hk) <- simple
n_hkLF <- melt(as.matrix(n_hk))
colnames(n_hkLF) <- c("Gene", "Sample", "Count")
ggplot(n_hkLF, aes(x=Gene, y=Count, color=Sample)) + geom_boxplot()
n_c <- n_counts[c(g_DEGorder, oarray, housekeep), ]
colnames(n_c) <- simple
n_DEGLF <- melt(as.matrix(n_c))
colnames(n_DEGLF) <- c("Gene", "Sample", "Count")
ggplot(n_DEGLF, aes(x=Gene, y=Count, color=Sample)) + geom_boxplot()

# Allomap ----
g_allomap <- c("Itgam", "Flt3", "Il1r2", "G6b", "Pf4", "Wdr40a", "Mir", "Arhu", 
               "Pdcd1", "Itga4", "Sema7a")
g_allomap <- c("Itgam", "Flt3", "Il1r2", "G6b", "Pdcd1", "Itga4", "Sema7a")
heatmap.2(t(t(n_order[g_allomap, ])), scale="row", col = coul, key=T, 
          xlab="Mouse", ylab="DESEQ2 Gene Expression", margins = c(6, 7), 
          Colv = F,
          key.title=NA, key.ylab=NA, trace="none", tracecol=NA)

# Database Genesets ----
# Need to run the T tests or at least filter the genes first!
# https://pubchem-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pathway/Reactome:R-MMU-8953897
# https://pubchem-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pathway/WikiPathways:WP88
# https://pubchem-ncbi-nlm-nih-gov.proxy.lib.umich.edu/pathway/Reactome:R-MMU-1280215
# Reactome Immunology Geneset
immun_set <- read.table(".\\Inputs\\2022 Tx Ideas\\ImmunologyPathwayID_12934_pcget_pathway_gene.txt", 
                        sep ="", header=T) # Raw Data
n_immun <- na.omit(n_addfilt[immun_set[,2],])
genes_immun <- rownames(n_immun)
coul <- colorRampPalette(brewer.pal(11, "RdBu"))(25)
colSide <- grouped
for (x in 1:sample_num) {
  if (grouped[x] == "Healthy") {
    colSide[x] <- "blue"
  } else if (grouped[x] == "Early_Rejecting") {
    colSide[x] <- "deeppink2"
  } else {
    colSide[x] <- "darkred" # Dark Red is Late rejecting
  }  }
heatmap.2(as.matrix(n_immun), scale = "row", col = coul, key = T, xlab = "Scaffold", 
          ylab="Reactome Immunity Geneset Expression", margins = c(7, 7),
          # density.info="none",
          key.title=NA, key.ylab=NA, trace="none", keysize=1.0, ColSideColors=colSide, 
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_immun), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=mouse, var.scale=1, circle=T) + 
  ggtitle("Reactome Immunity Geneset") + theme_classic()
# Reactome Adaptive Immunology Geneset
adapt_set <- read.table(".\\Inputs\\2022 Tx Ideas\\AdaptivePathwayID_12934_pcget_pathway_gene.txt", 
                        sep ="", header=T) # Raw Data
n_adapt <- na.omit(n_addfilt[adapt_set[,2],])
genes_adapt <- rownames(n_adapt)
heatmap.2(as.matrix(n_adapt), scale = "row", col = coul, key = T, xlab = "Scaffold", 
          ylab="Reactome Adaptive Immunity Geneset Expression", margins = c(7, 7),
          # density.info="none",
          key.title=NA, key.ylab=NA, trace="none", keysize=1.0, ColSideColors=colSide, 
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_adapt), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("Reactome Adaptive Immunity Geneset") + theme_classic()
# Reactome Innate Immunology Geneset
innate_set <- read.table(".\\Inputs\\2022 Tx Ideas\\InnatePathwayID_12934_pcget_pathway_gene.txt", 
                         sep ="", header=T) # Raw Data
n_innate <- na.omit(n_addfilt[innate_set[,2],])
genes_innate <- rownames(n_innate)
heatmap.2(as.matrix(n_innate), scale = "row", col = coul, key = T, xlab = "Scaffold", 
          ylab="Reactome Complement Cascade Geneset Expression", margins = c(7, 8),
          # density.info="none",
          key.title=NA, key.ylab=NA, keysize=1.0, ColSideColors=colSide, trace="none",
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_innate), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("Reactome Innate Immunity Geneset") + theme_classic()
# Landscape of innate immune system transcriptome and acute T cell-mediated 
# rejection of human kidney allografts
# Reactome TLR Cascade
tlr_set <- read.table(".\\Inputs\\2022 Tx Ideas\\TLRCascadePathwayID_12967_pcget_pathway_gene.txt", 
                      sep ="", header=T) # Raw Data
n_tlr <- na.omit(n_addfilt[tlr_set[,2],])
genes_tlr <- rownames(n_tlr)
heatmap.2(as.matrix(n_tlr), scale = "row", col = coul, key = T, xlab = "Scaffold", 
          ylab="Reactome Complement Cascade Geneset Expression", margins = c(7, 8),
          # density.info="none",
          key.title=NA, key.ylab=NA, keysize=1.0, ColSideColors=colSide, trace="none",
          Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_tlr), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("Reactome TLR Cascade Geneset") + theme_classic()
# Reactome PRR Geneset
# mRNAs encoding transmembrane and cytosolic PRRs
# mRNAs encoding caspases
# mRNAs encoding cytokines, chemokines, IFNs, and TNF
# mRNAs encoding secreted PRRs and the complement system proteins
comp_set <- read.table(".\\Inputs\\2022 Tx Ideas\\ComplementPathwayID_13006_pcget_pathway_gene.txt", 
                       sep ="", header=T) # Raw Data
n_comp <- na.omit(n_addfilt[comp_set[,2],])
genes_comp <- rownames(n_comp)
heatmap.2(as.matrix(n_comp), scale = "row", col = coul, key = T, xlab = "Scaffold", 
          ylab="Reactome Complement Cascade Geneset Expression", margins = c(7, 8),
          # density.info="none",
          key.title=NA, key.ylab=NA, keysize=1.0, ColSideColors=colSide, trace="none",
          # Colv = "Na",
          dendrogram = "none")
ggbiplot(prcomp(t(n_comp), scale.=T), ellipse=T, groups=simple, var.axes=F, 
         labels=samples, var.scale=1, circle=T) + 
  ggtitle("Reactome Complement Cascade Geneset") + theme_classic()
# mRNAs encoding DAMPs
# mRNAs encoding DNA damage sensors

# xCell ----
fac_DEseq <- grouped
coldata <- data.frame(cbind(fac_DEseq))
row.names(coldata) <- samples
dds <- DESeqDataSetFromMatrix(countData=r_addfilt, colData=coldata, design = ~ fac_DEseq)
paste(nrow(dds), " genes input into DESeq2 for GSEA", sep="")
dds <- DESeq(dds)
human_genes1 <- convertMouseGeneList(as.array(rownames(get_counts(dds))))
human_genes1 <- human_genes1[!duplicated(human_genes1[1]),]
human_genes2 <- as.data.frame(human_genes1[,2])
row.names(human_genes2) <- human_genes1[,1]
colnames(human_genes2) <- "human_genes"
combined_list <- merge(get_counts(dds), human_genes2, by = "row.names", all.x = T)
nulllist <- combined_list[is.na(combined_list$human_genes),]
nulllist$human_genes <- nulllist$Row.names
notnulllist <- combined_list[!is.na(combined_list$human_genes),]
combined_list2 <- rbind(nulllist, notnulllist)
write.csv(combined_list2, file= "2022_02_24_STx_t-test_for_xCell.csv")
# Input into xCell  https://xcell.ucsf.edu/
# View Results
xC_res <- read.table(".\\Inputs\\2022 Tx Ideas\\xCell_2022_02_24_STx_t-test_for_xCell_xCell_0551022422.txt", sep ="", header=T)
xC_scr <- read.table(".\\Inputs\\2022 Tx Ideas\\xCell_2022_02_24_STx_t-test_for_xCell_xCell_0551022422_RAW.txt", sep ="", header=T)
xC_pval <- read.table(".\\Inputs\\2022 Tx Ideas\\xCell_2022_02_24_STx_t-test_for_xCell_xCell_0551022422.pvals.txt", sep ="", header=T)