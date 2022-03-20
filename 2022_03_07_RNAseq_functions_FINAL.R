# Sample Differences ----
# To briefly visualize which samples are most alike and most different.
sample.diff <- function(counts_df, sample_names) {
  mat_dist <- as.matrix(dist(t(counts_df))) # Euclidean dist between samples
  mat_dist <- mat_dist/max(mat_dist)# Normalized by the max distance
  cim(mat_dist, symkey=F, row.names=sample_names, col.names=sample_names)
}

# MA Plots ----
plot.topMA <- function(counts_df) {
  sum_MA <- matrix(NA, nrow=sample_num, ncol=sample_num) # Empty matrix for loop
  pb <- txtProgressBar(min=0,max=(sample_num^2)-ncol(combn(sample_num, 2)), 
                       style=3) # Progress bar
  k <- 0 # Progress bar counter, starts at 0
  # WARNING: this nested loop is HELLA slow! (approx 30 sec per sample)
  num_samples <- ncol(counts_df)
  for (i in 1:num_samples) {
    for (j in i:num_samples) {
      k <- k + 1 # Progress bar counter
      setTxtProgressBar(pb, k) # Update progress bar
      if (i!=j){
        M <- counts_df[,i] - counts_df[,j]
        A <- (counts_df[,i] + counts_df[,j])/2
        df <- data.frame(A, M)
        sum_MA[j,i] <- sum(abs(predict(loess(M~A,df), df$A)))}}}
  # Which paired samples are most different and create MA plots of the top 9.
  num_plots <- 9 # In case you're bored enough to want more MA plots
  topMA <- which(sum_MA>=sort(sum_MA, decreasing=T)[num_plots], arr.ind=T)
  topMA <- topMA[order(sum_MA[topMA], decreasing=T),]
  plot_lst <- c() # Empty plot list for loop
  for (i in seq(1,num_plots)) {
    M <- counts_df[,topMA[i,1]] - counts_df[, topMA[i,2]]
    A <- (counts_df[,topMA[i,1]] + counts_df[, topMA[i,2]])/2
    df <- data.frame(A, M)
    plot_lst[[i]] <- ggplot(df, aes(x = A, y = M)) + 
      ylim(-6,6) + # ylim throws an error, but still works
      geom_hline(yintercept=0, color="blue3") + # Blue line through y=0
      geom_point(size = 0.8, alpha = 1/5) + 
      stat_smooth(se = F, method = "loess", color = "red3") + # Red smoothed line
      ggtitle(paste(samples[topMA[i,1]], "vs.", samples[topMA[i,2]]))}
  ggarrange(plotlist=plot_lst, common.legend=F, nrow=3, ncol=3)
}

# PCA ----
plot.pca <- function(counts_df, fac, ptitle="", ellipses=F) {
  pca_data <- prcomp(counts_df, center=T, scale=T) # Calculate PCA
  d <- data.frame(cbind(pca_data$rotation[,1],pca_data$rotation[,2]))
  g <- ggplot(d, aes(x=X1, y=X2, color=fac)) + 
    geom_point(size=3, shape = 19) + 
    labs(title=paste("PCA, ", ptitle, sep=""), x="PC1", y="PC2", color="Groups") + 
    theme_classic() + theme(legend.position = "right")
  if (ellipses==T) {
    g <- g + ggforce::geom_mark_ellipse(aes( color = fac), 
                                        show.legend=F)
  }
  g
} 

# Scree Plot ----
plot.scree <- function(pca_data, ptitle="") {
  pca_var <- ((pca_data$sdev^2)/(sum(pca_data$sdev^2)))*100 # Proportion of variance
  x <- barplot(pca_var, xlab=paste("Principal Component, 1 -", length(pca_data$sdev)),
               ylab="Proportion of variance (%)", names.arg=c(1:length(pca_data$sdev)),
               main=paste("Scree plot, ", ptitle, sep=""),
               ylim=c(0,102)) # ylim=102 shows text
  y <- round(pca_var, digits=2) # Round values
  text(x,y,labels=as.character(paste(y,"%", sep="")), pos=3, cex=0.6) # Add values
}

# Visualize Variation ----
plot.variation <- function(lst.hist, line_types, sample_names, ylim1=2000, ylim2=2000) {
  num_samples <- length(lst.hist)
  # Remove genes with no expression
  lst.hist_nz <- lapply(1:num_samples, function (x) lst.hist[[x]][lst.hist[[x]][,2]>0,])
  avg <- apply(simplify2array(lst.hist), 1:2, mean) # Average of the histograms
  avg_nz <- avg[avg[,2]>0,] # Remove zeros from average
  layout(matrix(c(1,1,2,3), ncol=2, byrow=F)) # All on the same plot
  par(mar = c(4,4,2,1), oma = c(1,1,1,1)) # Adjust margins
  plot(lst.hist_nz[[1]][,1],lst.hist_nz[[1]][,2], type="l", frame=T, pch=4, col="red", 
       xlab=expression(log[2](count + 1)), ylab="Frequency", 
       xlim=c(0,15), ylim=c(0,ylim1), oma=c(0,0,0,0), mar=c(0,0,0,0)) # Plot 1st line
  axis(1, seq(0,15,1)) # Specify axis ticks
  legend("topright", legend=sample_names, col=coul, lty=line_types, yjust=1)
  for(i in 2:num_samples){ # Add additional lines
    lines(lst.hist_nz[[i]][,1],lst.hist_nz[[i]][,2], pch=4, col=coul[i], type="l", 
          lty=line_types[i])}
  plot(lst.hist_nz[[1]][,1],lst.hist_nz[[1]][,2], type="l", frame=T, pch=4, col="red", # Zoom1
       xlab="", ylab="", xlim=c(0,3), ylim=c(0,ylim2), main="Zoom1", 
       oma=c(0,0,0,0), mar=c(0,0,0,0)) # xlim and ylim set MANUALLY
  for(i in 2:num_samples){ # Add additional lines
    lines(lst.hist_nz[[i]][,1],lst.hist_nz[[i]][,2], pch=4, col=coul[i], type="l", 
          lty=line_types[i])}
  plot(lst.hist_nz[[1]][,1],lst.hist_nz[[1]][,2], type="l", frame=T, pch=4, col="red", # Zoom2
       xlab="", ylab="", xlim=c(3,10), ylim=c(0,200), main="Zoom2", 
       oma=c(0,0,0,0), mar=c(0,0,0,0)) # xlim and ylim set MANUALLY
  for(i in 2:num_samples){ # Add additional lines
    lines(lst.hist_nz[[i]][,1],lst.hist_nz[[i]][,2], pch=4, col=coul[i], type="l", 
          lty=line_types[i])}
  mtext("Finding Filters, Pseudo Counts", side=3, line=-0.5, outer=T) # Add title
}

# Visualize Variation Short----
plot.variationshort <- function(lst.hist, line_types, sample_names) {
  num_samples <- length(lst.hist)
  # Remove genes with no expression
  lst.hist_nz <- lapply(1:num_samples, function (x) lst.hist[[x]][lst.hist[[x]][,2]>0,])
  avg <- apply(simplify2array(lst.hist), 1:2, mean) # Average of the histograms
  avg_nz <- avg[avg[,2]>0,] # Remove zeros from average
  plot(lst.hist_nz[[1]][,1],lst.hist_nz[[1]][,2], type="l", frame=T, pch=4, col="red", 
       xlab=expression(log[2](count + 1)), ylab="Frequency", 
       oma=c(0,0,0,0), mar=c(0,0,0,0)) # Plot 1st line
  axis(1, seq(0,15,1)) # Specify axis ticks
  legend("topright", legend=sample_names, col=coul, lty=line_types, yjust=1)
  for(i in 2:num_samples){ # Add additional lines
    lines(lst.hist_nz[[i]][,1],lst.hist_nz[[i]][,2], pch=4, col=coul[i], type="l", 
          lty=line_types[i])}
  mtext("Pseudo Filtered", side=3, line=-0.5, outer=T) # Add title
}

# Error from Mean ----
plot.error <- function(lst.hist, line_types) {
  num_samples <- length(lst.hist)
  avg <- apply(simplify2array(lst.hist), 1:2, mean) # Average of the histograms
  avg_nz <- avg[avg[,2]>0,] # Remove zeros from average
  layout(matrix(c(1,1,2,2), ncol=2, byrow=F)) # All on the same plot
  par(mar = c(4,4,2,2), oma = c(1,1,1,1)) # Adjust margins to better use space
  sum <- 0 # Reset the sum
  for(i in 1:num_samples){ # Sum bins for each count
    sum <- sum + abs(lst.hist[[i]][,2]-avg[,2])}
  sum <- cbind(lst.hist[[1]][,1],sum)
  sum_nz <- sum[sum[,2]>0,] # Remove zeros
  plot(sum_nz, xlab=expression(log[2](count + 1)), ylab="Sum of Differences",
       xlim=c(0,15), oma=c(0,0,0,0), mar=c(0,0,0,0), 
       main="Sum of Differences from the Mean", cex=0.75)
  axis(1, seq(0,15,1)) # Specify axis ticks
  lines(sum_nz, type="l") # Connect points
  sum_nz_norm <- cbind(sum_nz[,1], (sum_nz[,2]/avg_nz[,2])) # Normalized error
  sum_nz_norm <- na.omit(sum_nz_norm)
  plot(sum_nz_norm, xlab=expression(log[2](count + 1)), xlim=c(0,15), 
       oma=c(0,0,0,0), mar=c(0,0,0,0), ylab="Sum of Differences/Avg",
       main="Sum of Differences from the Mean/Avg", cex=0.75)
  lines(sum_nz_norm, type="l")
  axis(1, seq(0,15,1)) # Specify axis ticks
  mtext("Finding Filters, Pseudo Counts", side=3, line=-0.5, outer=T) # Add title
}


# Test Filters ----
filters.genes <- function( counts_df, low_filters, high_filters, DESeq_factor) {
  
  sample_names <- colnames(counts_df)
  num_lowfilt <- length(low_filters)
  num_highfilt <- length(high_filters)
  
  g_filtered <- lapply(1:num_highfilt, function (x) {
    keep <- rowSums(counts_df) < (high_filters[x]*sample_num) # High count filter
    counts_filt <- counts_df[keep,]
    lapply(1:num_lowfilt, function (y) {
      keep <- rowSums(counts_filt) > (low_filters[y]*sample_num) # Low count filter
      counts_filt <- counts_filt[keep,]
      g_filt <- rownames(counts_filt) # Vector of filtered genes
    })
  })
  
  names(g_filtered) <- paste(high_filters, " high filter", sep="")
  
  for (i in 1:length(g_filtered)) {
    names(g_filtered[[i]]) <- paste(" ", low_filters, " low filter", sep="")
  }
  g_filtered <- unlist(g_filtered,recursive=FALSE)

  # Filtered raw counts used for DESeq2
  r_filtered <- lapply(1:length(g_filtered), function (x) r_counts[g_filtered[[x]],])
  names(r_filtered) <- names(g_filtered)
  
  # Filtered cpm used for plotting
  n_filtered <- lapply(1:length(g_filtered), function (x) n_counts[g_filtered[[x]],])
  names(n_filtered) <- names(g_filtered)
  
  # Filtered pseudo norm counts used for ??????????????????????
  ps_filtered <- lapply(1:length(g_filtered), function (x) ps_counts[g_filtered[[x]],])
  names(ps_filtered) <- names(g_filtered)
  
  return(list("FilteredGenes"=g_filtered, "RawFiltered"=r_filtered, 
              "CPMFiltered"=n_filtered, "PseudoFiltered"=ps_filtered))
}

# Find DEGS from filters.genes ----
filters.DEGs <- function( filters.genes_list, DESeq_factor) {
  coldata <- data.frame(cbind(DESeq_factor)) # DEG Analysis
  raw_filtered <- filters.genes_list$RawFiltered
  row.names(coldata) <- colnames(raw_filtered[[1]])
  pb <- txtProgressBar(min = 0, max = length(raw_filtered), style = 3)
  DE_genes <- lapply(1:length(raw_filtered), function(x) {
    setTxtProgressBar(pb, x)
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData=raw_filtered[[x]], 
                                  colData=coldata, design = ~ DESeq_factor))
    dds <- suppressMessages(DESeq(dds))
    res <- results(dds) # Table of log2 fold changes, p-values, & p-adj values
    res_ordered <- res[order(abs(res$log2FoldChange)),] # Order results by LFC
    dds_genes <- rownames(subset(res_ordered, padj < 0.1)) # Subset by p-adj < 0.1
  })
  names(DE_genes) <- names(raw_filtered)
  return(list("DEGs"=DE_genes))
}

# SVD ----
filters.SVDs <- function( filters.genes_list, DESeq_factor) { # Function needs to be adjusted for each case!
  
  n_DEG <- filters.genes_list$CPMFiltered
  sample_names <- colnames(filters.genes_list$CPMFiltered[[1]])
  sample_number <- length(sample_names)
  healthy_samples <- DESeq_factor == "Healthy"
  
  SVD_euc <- lapply(1:length(n_DEG), function (x) {
    centered <- n_DEG[[x]] - rowMeans(n_DEG[[x]]) # First, center data on genes
    svd1 <- svd(t(centered))
    Z <- t(centered) %*% svd1$v
    centroid_healthy <- colMeans(rbind(svd1$u[healthy_samples,]),dims=1)
    euclid <- rbind(centroid_healthy, svd1$u)
    euclid_dist <- rdist(euclid[,1:3], metric = "euclidean", p = 2)
    euclid_dist <- euclid_dist[1:sample_number]
    euc <- cbind.data.frame(sample_names, euclid_dist)
    rownames(euc) <- c(sample_names)
    euc_SVD <- cbind(euc, DESeq_factor)
  })
  names(SVD_euc) <- names(filters.genes_list$CPMFiltered)
  return(list("SVD"=SVD_euc))
}

# RF ----
filters.RF <- function( filters.genes_list, DESeq_factor, clinical_info, study_1_term) {
  n_DEG <- filters.genes_list$CPMFiltered
  sample_names <- colnames(filters.genes_list$CPMFiltered[[1]])
  sample_number <- length(sample_names)
  healthy_samples <- DESeq_factor == "Healthy"
  pb <- txtProgressBar(min = 0, max = length(n_DEG), style = 3)
  RF_euc <- lapply(1:length(n_DEG), function(x) {
    setTxtProgressBar(pb, x)
    predictor_data <- t(n_DEG[[x]]) # Load Data as Predictor Variable
    target <- clin_info[,"state"] # Set Target Variable
    target[target==0] <- "Healthy"
    target[target==1] <- study_1_term
    target <- as.factor(target)
    tmp <- as.vector(table(target)) # Run RF Algorithm
    num_classes <- length(tmp)
    min_size <- tmp[order(tmp,decreasing=F)[1]]
    sampsizes <- rep(min_size,num_classes)
    rf_output <- randomForest(x=predictor_data, y=target, importance=T, ntree=10001, 
                              proximity=T, sampsize=sampsizes, na.action=na.omit)
    plot_MDS <- MDSplot(rf_output, target, k=2) # MDS plot: Separation of Classes
    plot_MDS <- as.data.frame(plot_MDS$points)
    centroid_healthy <- colMeans(rbind(plot_MDS[healthy_samples,]),dims=1)
    euclid <- rbind(centroid_healthy, plot_MDS) # Euclidean distances
    euclid_dist <- rdist(euclid[,1:2], metric = "euclidean", p = 2)
    euclid_dist <- euclid_dist[1:sample_number]
    euc <- cbind.data.frame(sample_names, euclid_dist)
    rownames(euc) <- c(sample_names)
    euc_RF <- cbind(euc, DESeq_factor)
  })
  names(RF_euc) <- names(filters.genes_list$CPMFiltered)
  return(list("RF"=RF_euc))
}

# Scores ----
filters.scores <- function( SVD_list, RF_list) {
  sample_names <- SVD_list$SVD[[1]]$sample_names
  euc_SVD <- SVD_list$SVD
  euc_RF <- RF_list$RF
  scores <- lapply(1:length(euc_SVD), function(x) {
    scoring <- cbind.data.frame(euc_SVD[[x]]$euclid_dist, euc_RF[[x]]$euclid_dist)
    colnames(scoring) <- c("SVD", "RF")
    rownames(scoring) <- sample_names
    scoring
  })
  names(scores) <- names(euc_SVD)
  return(list("Scores"=scores))
}

# Plot Scores ----
filters.scoreplots <- function( scores_list, color_factor, shape_factor, label_factor) {
  filter_sets <- names(scores_list$Scores)
  plot_lst <- lapply(1:length(scores_list$Scores), function (x) {
    ggplot(scores_list$Scores[[x]], aes(x = SVD,y = RF, color = color_factor)) + 
      geom_point(aes(shape = shape_factor), size = 1) + 
      geom_text(aes(label = label_factor), nudge_x = 0.03, nudge_y = -0.01, size = 3) + 
      xlab("SVD (arb. units)") + ylab("RF (probability)") + 
      ggtitle(filter_sets[x]) + 
      ggforce::geom_mark_ellipse(aes(fill = color_factor, color = color_factor)) + 
      coord_equal()
  })
  names(plot_lst) <- filter_sets
  return(list("ScorePlots"=plot_lst))
}

# Table of Sig P vals ----
sig_pvaltable <- function(entries, entries_names) {
  # Creates a table of the number of genes with p values < 0.05 for EVERY pair
  # of the sample data groups. Warning: This function takes a LONG TIME to run 
  # because it generates a huge number of t tests.
  num_genes <- nrow(data.frame(entries[1])) # Number of genes for each t test
  num_entries <- length(entries) # Number of groups for use in t tests
  sig_pvals <- data.frame(matrix(NA, nrow = num_entries, ncol = num_entries))
  # Empty table to fill
  rownames(sig_pvals) <- entries_names
  colnames(sig_pvals) <- entries_names
  pb <- txtProgressBar(min = 0, max = (num_entries^2-ncol(combn(num_entries, 2))), 
                       style = 3) # Progress bar
  l <- 0
  for (i in 1:num_entries){
    for (j in i:num_entries){
      if (i != j) { # Don't test identical groups
        pval <- vector("double", num_genes) # Pre-allocate results
        cond1 <- data.frame(entries[i])
        cond2 <- data.frame(entries[j])
        try(for (k in seq(num_genes)){
          pval[k] <- t.test(cond1[k,], cond2[k,])$p.value}) # T test
        pval <- data.frame(pval) # As data frame
        pval[is.na(pval)] <- 0 # Replace NaN's with zeros
        sig_pvals[i,j] <- sum(pval > 0 & pval <= 0.05) # Significant p values
      } else {sig_pvals[i,j] <- NA} # No point in a t test of identical groups
      l <- l + 1 # Progress bar
      setTxtProgressBar(pb, l)}} # Progress bar
  sig_pvals} # Output table showing number of sig. p values for each pair

# pvalgenes Function ----
pvalgenes <- function(group1, group2, num_breaks=100, plot=F) {
  # Plots the distribution of p values for a pair of sample groups. P values 
  # < 0.05 are colored green.
  pval <- sapply(1:nrow(group1), function(x) { # Try function sidesteps errors
    try(t.test(group1[x,], group2[x,])$p.value)})
  pval <- data.frame(pval) # As data frame
  rownames(pval) <- rownames(group1)
  if (plot == T) { # For plotting the distribution
    p_hist <- hist(pval[pval>0], breaks=num_breaks, plot=F) # Store histogram
    coul <- ifelse(p_hist$breaks<=0.05, "light green", "light gray") # Colors
    p_title <- paste(deparse(substitute(group1)), " & ", deparse(substitute(
      group2)), ", ", sum(pval > 0 & pval <= 0.05), " genes", sep="") # Title
    plot(p_hist, col=coul, border=F, main=p_title, xlab="p value", ylab="Freq")}
  pval <- pval[pval != 0 & pval <= 0.05, , drop = F]
  rownames(pval)}

# plot_pca Function ----
# This function is still a WORK IN PROGRESS, but it does work
plot_pca <- function(gene_counts) {
  coul <- c(1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,3,3,3)
  pca <- prcomp(gene_counts, center=T, scale=T)
  plot(pca$rotation[,1],pca$rotation[,2], xlab = "PC1", ylab = "PC2",
       col = coul, pch = 19, cex = 1, 
       main=paste(deparse(substitute(gene_counts))))
  text(pca$rotation[,1],pca$rotation[,2], samples, cex=0.7, pos=4)
  dataEllipse(pca$rotation[c(1:9,13:15),1], pca$rotation[c(1:9,13:15),2], 
              levels=c(0.7), center.pch=F, draw=T, add=T, segments=51, 
              robust=F, plot.points=F, col=1, pch=1, lwd=1, lty=2)
  dataEllipse(pca$rotation[c(10:12),1], pca$rotation[c(10:12),2], 
              levels=c(0.7), center.pch=F, draw=T, add=T, segments=51, 
              robust=F, plot.points=F, col=2, pch=1, lwd=1, lty=2)
  dataEllipse(pca$rotation[c(16:18),1], pca$rotation[c(16:18),2], 
              levels=c(0.7), center.pch=F, draw=T, add=T, segments=51, 
              robust=F, plot.points=F, col=3, pch=1, lwd=1, lty=2)
  summary(pca)}

# Yining's GSEA Functions ----
get_counts <- function(x){
  normcounts <- counts(x, normalized = TRUE)
  return (normcounts)
}

# write results to a table
write_results <- function(tableX, table_nameX){
  write.csv(as.data.frame(tableX), file = table_nameX)
}

# get results table from a written file
get_results <- function(table_nameX){
  results_table <- read.table(table_nameX, sep =",", header = FALSE)
  return(results_table)
}

convertMouseGeneList <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <- getLDS(attributesL = c("hgnc_symbol"),
                    filters = "mgi_symbol",
                    values = x,
                    mart = mouse,
                    attributes = c("mgi_symbol"),
                    martL = human,
                    uniqueRows=FALSE)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  human_genes <- genesV2
  return(human_genes)
}

# Elastic Net ----
multinom3_EN <- function(x, y, num_folds){
  alphaList <- seq(1.0, 0.0, -0.05)
  pb <- txtProgressBar(min = 0, max = length(alphaList), style = 3)
  ENgenes <- lapply(1:length(alphaList), function(i) {
    setTxtProgressBar(pb, i)
    cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alphaList[i], 
                                        nfolds = num_folds, family = "multinomial"))
    best.lambda <- cvfit$lambda.min
    fit <- suppressMessages(glmnet(x, y, alpha = alphaList[i], lambda = best.lambda, 
                  family = "multinomial"))
    Coef<-coef(fit, s = best.lambda)
    Index<-c(Coef[["0"]]@i[-1], Coef[["1"]]@i[-1], Coef[["2"]]@i[-1])
  })
  
  names(ENgenes) <- alphaList
  return(list("ENgenes"=ENgenes))
 
}

binom_EN <- function(x, y, num_folds){
  alphaList <- seq(1.0, 0.0, -0.05)
  pb <- txtProgressBar(min = 0, max = length(alphaList), style = 3)
  ENgenes <- lapply(1:length(alphaList), function(i) {
    setTxtProgressBar(pb, i)
    cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alphaList[i], 
                                        nfolds = num_folds, family = "multinomial"))
    best.lambda <- cvfit$lambda.min
    fit <- suppressMessages(glmnet(x, y, alpha = alphaList[i], lambda = best.lambda, 
                                   family = "multinomial"))
    Coef<-coef(fit, s = best.lambda)
    Index<-c(Coef[["0"]]@i[-1], Coef[["1"]]@i[-1])
  })
  
  names(ENgenes) <- alphaList
  return(list("ENgenes"=ENgenes))
  
}

# Old code for EN with a single alpha value, alpha=1.0 LASSO
# cvfit <- cv.glmnet(x, y, alpha=1.0, nfolds=100,
#                    family = "multinomial") # simple=binomial, grouped=multinomial
# plot(cvfit)
# cvfit$lambda.min
# best.lambda <- cvfit$lambda.min
# fit <- glmnet(x, y, 
#               alpha=1.0, 
#               lambda = best.lambda, 
#               family = "multinomial")
# plot(fit, xvar = "lambda", label = T, type.coef = "2norm")
# print(fit)
# coef(fit, s = best.lambda)
# Coef<-coef(fit, s = best.lambda)
# Coef[["0"]]@i[-1]
# Coef[["1"]]@i[-1]
# Coef[["2"]]@i[-1]
# Index<-c(Coef[["0"]]@i[-1],
#          Coef[["1"]]@i[-1],
#          Coef[["2"]]@i[-1]
#          )
# Index
# genes_ENL <- rownames(r_ttest)[Index]
# genes_ENL

# GSEA Analyses ----
plot_GSEA <- function(folder, name_dataset, num_genesets){
  # Load NES Scores
  namefile1 <- paste(folder, "\\RvH.tsv", sep = "")
  namefile2 <- paste(folder, "\\EvH.tsv", sep = "")
  namefile3 <- paste(folder, "\\LvH.tsv", sep = "")
  namefile4 <- paste(folder, "\\EvL.tsv", sep = "")
  nes_table1 <- read.table(namefile1, sep = "\t", header = T, fill = T)
  nes_table1 <- nes_table1[, c(-2,-3,-12)]
  nes_table2 <- read.table(namefile2, sep = "\t", header = T, fill = T)
  nes_table2 <- nes_table2[, c(-2,-3,-12)]
  nes_table3 <- read.table(namefile3, sep = "\t", header = T, fill = T)
  nes_table3 <- nes_table3[, c(-2,-3,-12)]
  nes_table4 <- read.table(namefile4, sep = "\t", header = T, fill = T)
  nes_table4 <- nes_table4[, c(-2,-3,-12)]
  # DFs of rows=gene sets
  gs1 <- nes_table1[,1] # Extract first 40 characters from gene set name
  gs2 <- nes_table2[,1]
  gs3 <- nes_table3[,1]
  gs4 <- nes_table4[,1]
  nes1 <- nes_table1[,-c(2,3,7,8,9)]
  nes2 <- nes_table2[,-c(2,3,7,8,9)]
  nes3 <- nes_table3[,-c(2,3,7,8,9)]
  nes4 <- nes_table4[,-c(2,3,7,8,9)]
  nes1 <- nes1[1:num_genesets,] # Top NES gene sets
  nes2 <- nes2[1:num_genesets,] # Top NES gene sets
  nes3 <- nes3[1:num_genesets,] # Top NES gene sets
  nes4 <- nes4[1:num_genesets,] # Top NES gene sets
  genesets1 <- nes1[,1] # Extract first 40 characters from gene set name
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
  tab1 <- nes_table1[genesets, -c(2,3,7,8,9)]
  tab2 <- nes_table2[genesets, -c(2,3,7,8,9)]
  tab3 <- nes_table3[genesets, -c(2,3,7,8,9)]
  tab4 <- nes_table4[genesets, -c(2,3,7,8,9)]
  nes_top <- cbind.data.frame(tab1, tab2[,2:4], tab3[,2:4], tab4[,2:4])
  colnames(nes_top) <- c("GENESET", "Rejecting v Healthy", "p1", "q1", 
                         "Early Rejecting v Healthy", "p2", "q2", 
                         "Late Rejecting v Healthy", "p3", "q3", 
                         "Early Rejecting v Late", "p4", "q4")
  write.table(genesets, file = paste(name_dataset, ".csv", sep = ""), sep="")
  genesets <- substr(genesets, 1, 60) # Extract first 60 characters from gene set name
  nes_top[,1] <- genesets
  rownames(nes_top) <- genesets
  nes_top[is.na(nes_top)] = 0
  nes_p <- nes_top[,-c(2,4,5,7,8,10,11,13)]
  nes_pLF <- melt(nes_p, id = c("GENESET"))
  nes_top <- nes_top[,-c(3,4,6,7,9,10,12,13)]
  nes_topLF <- melt(nes_top, id = c("GENESET"))
  nes_topLF$variable <- factor(nes_topLF$variable,levels=unique(nes_topLF$variable))
  nes_topLF <- cbind.data.frame(nes_topLF, nes_pLF[,3])
  colnames(nes_topLF) <- c("GENESET", "variable", "NES", "pval")
  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57")
  xx <- ggplot(nes_topLF, aes(x = GENESET, y = variable, label=round(NES, digits=2))) + 
    geom_point(aes(size = NES, fill = -log10(pval)), alpha = 0.75, shape = 21) + 
    geom_text(size=3, nudge_y = 0.25, color = ifelse(nes_topLF$NES < 1.50, 1, 2)) + 
    scale_color_manual(values = c("deepred", "black")) + 
    labs( x= "Geneset", y = "Comparison", size = "NES", fill = "-log10 pval")  + 
    # scale_fill_manual(values = colours, guide = FALSE) + 
    # theme_light() +
    theme(axis.text.x  = element_text(angle=30, hjust=0.95,vjust=1.05)) + coord_flip() + 
    theme(panel.grid.major.y = element_blank(), panel.border = element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.ticks.x = element_blank(),
          # axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  xx

}

# Old Code Ideas
# ssizess <- (-log10(nes_topLF$pval)+2)^2
# ggplot(nes_topLF, aes(x=GENESET, y=NES, label=round(NES ,digits=2), color=NES)) + 
#   geom_segment(aes(y = 1.5, x = GENESET, yend = NES, xend = GENESET), color = "black") +
#   geom_point(stat='identity', fill="black", size=ssizess) + 
#   geom_text(color="black", size=3, nudge_y = nes_top$NES/50) +
#   labs(title="Top Gene Sets by NES", 
#        subtitle=paste("Skin Transplant Training Cohort ", dataset)) + 
#   scale_colour_gradient(low="yellow", high="red") + 
#   ylim(1.5, 2.25) +
#   theme(axis.title.y = element_text(size=10),
#         axis.text.y  = element_text(angle=0, size=10)) + 
#   theme_light() + coord_flip() +
#   theme(panel.grid.major.y = element_blank(), panel.border = element_blank(),
#         # axis.ticks.y = element_blank(),
#         # axis.ticks.x = element_blank(),
#         # axis.title.x = element_blank(),
#         # axis.title.y = element_blank(),
#         # axis.text.x = element_blank(),
#         # panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank())


# GGBIPLOT ----
# 
#  ggbiplot.r
#  
#  Copyright 2011 Vincent Q. Vu.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 

#' Biplot for Principal Components using ggplot2
#'
#' @param pcobj           an object returned by prcomp() or princomp()
#' @param choices         which PCs to plot
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group?
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param var.axes        draw arrows for the variables?
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#'
#' @return                a ggplot2 plot
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))
#'
ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, 
                     groups = NULL, ellipse = FALSE, ellipse.prob = 0.68, 
                     labels = NULL, labels.size = 3, alpha = 1, 
                     var.axes = TRUE, 
                     circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, 
                     varname.abbrev = FALSE, ...)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  stopifnot(length(choices) == 2)
  
  # Recover the SVD
  if(inherits(pcobj, 'prcomp')){
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$rotation
  } else if(inherits(pcobj, 'princomp')) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- pcobj$loadings
  } else if(inherits(pcobj, 'PCA')) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
  } else if(inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  } else {
    stop('Expected a object of class prcomp, princomp, PCA, or lda')
  }
  
  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))
  
  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])
  
  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)
  
  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  
  # Scale the radius of the correlation circle so that it corresponds to 
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  
  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }
  
  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs, 
                       sprintf('(%0.1f%% explained var.)', 
                               100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }
  
  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()
  
  if(var.axes) {
    # Draw circle
    if(circle) 
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + geom_path(data = circle, color = muted('white'), 
                         size = 1/2, alpha = 1/3)
    }
    
    # Draw directions
    g <- g +
      geom_segment(data = df.v,
                   aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = arrow(length = unit(1/2, 'picas')), 
                   color = muted('red'))
  }
  
  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    } else {
      g <- g + geom_text(aes(label = labels), size = labels.size)      
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    } else {
      g <- g + geom_point(alpha = alpha)      
    }
  }
  
  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    
    ell <- ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'), 
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  
  # Label the variable axes
  if(var.axes) {
    g <- g + 
      geom_text(data = df.v, 
                aes(label = varname, x = xvar, y = yvar, 
                    angle = angle, hjust = hjust), 
                color = 'darkred', size = varname.size)
  }
  # Change the name of the legend for groups
  # if(!is.null(groups)) {
  #   g <- g + scale_color_brewer(name = deparse(substitute(groups)), 
  #                               palette = 'Dark2')
  # }
  
  # TODO: Add a second set of axes
  
  return(g)
}

# 
#  ggscreeplot.r
#
#  Copyright 2011 Vincent Q. Vu.
# 
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 

#' Screeplot for Principal Components
#'
#' @param pcobj          an object returned by prcomp() or princomp()
#' @param type           the type of scree plot.  'pev' corresponds proportion of explained variance, i.e. the eigenvalues divided by the trace. 'cev' corresponds to the cumulative proportion of explained variance, i.e. the partial sum of the first k eigenvalues divided by the trace.
#' @export
#' @examples
#'   data(wine)
#'   wine.pca <- prcomp(wine, scale. = TRUE)
#'   print(ggscreeplot(wine.pca))
#'
ggscreeplot <- function(pcobj, type = c('pev', 'cev')) 
{
  type <- match.arg(type)
  d <- pcobj$sdev^2
  yvar <- switch(type, 
                 pev = d / sum(d), 
                 cev = cumsum(d) / sum(d))
  
  yvar.lab <- switch(type,
                     pev = 'proportion of explained variance',
                     cev = 'cumulative proportion of explained variance')
  
  df <- data.frame(PC = 1:length(d), yvar = yvar)
  
  ggplot(data = df, aes(x = PC, y = yvar)) + 
    xlab('principal component number') + ylab(yvar.lab) +
    geom_point() + geom_path()
}
