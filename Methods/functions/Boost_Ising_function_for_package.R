# load packages

library(mclust)
library(lattice)
library(edgeR)
library(ggplot2)

source("functions/functions_potts.R");
Rcpp::sourceCpp('functions/functions_potts_omega.cpp');

transform_data <- function(count, loc){
  array_y <- array(0, dim = c((max(loc[,1])-min(loc[,1])+1),(max(loc[,2])-min(loc[,2])+1),ncol(count)))
  loc1 <- min(loc[,1])-1
  loc2 <- min(loc[,2])-1
  for (i in 1:ncol(count)){
    for (j in 1:nrow(loc)){
      array_y[loc[j,1]-loc1,loc[j,2]-loc2, i] <- count[j,i]
    }
  }
  return(array_y)
}


# filter data
filter.st<- function(count, sample_info, min_total = 10,min_percentage = 0.1){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  if(sum(rowSums(count) < min_total) == 0){
    sample_f <- sample_info
    count_f <- count[,-(which(colSums(count == 0) > (1-min_percentage)*sample_num))]
  }
  else{
  sample_f <- sample_info[-which(rowSums(count)<min_total),]
  count_f <- count[-which(rowSums(count)<min_total),-(which(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num))]
  }
  return(list(sample_f,count_f))
}


# plot ST data
plot.st <- function(count, loc, gene_name = NULL, main = NULL) {
  if (is.null(colnames(count))){stop('Please provide gene names as column names of count matrix')}
  if (is.null(gene_name) == FALSE){
    if (gene_name %in% colnames(count)){
      count <- count[, gene_name]
    } else {stop('Please enter a valid gene name')}
  }
  else{
    gene_name <- colnames(count)[which(colSums(count) == max(colSums(count)))]
    count <- count[, gene_name]
  }
  if (is.null(main)){main <- gene_name}
  data <- data.frame(expr = as.vector(count), x = loc[, 1], y = loc[, 2]);
  ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = 4) + 
    coord_fixed(ratio = 1) + scale_color_distiller(palette = "Spectral") + 
    theme_classic() + labs(color = "", title = main) + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5),
          legend.position="right")
}

# plot ST data
plot.st.binary <- function(count, loc, main = NULL) {
  data <- data.frame(expr = as.factor(as.vector(count)), x = loc[, 1], y = loc[, 2]);
  ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = 4) + 
    coord_fixed(ratio = 1) + scale_color_manual(values = c("1" = "blue", '2' = 'red'))  + 
    theme_classic() + labs(color = "", title = main) + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5),
          legend.position="right")
}

# get size factor estimate
get.size.factor <- function(count, norm_method = 'tss'){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors
    count_q75 <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    count_N <- rowSums(count)/nrow(count)
    raw_s_factors <- count_q75/count_N
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    count_t <- t(count)
    raw_s_factors <- calcNormFactors(count_t,method = "TMM")
  }
  else{
    stop("Please choose a valid normalization method")
  }
  return(raw_s_factors)
}

# st data normalization
normalize.st <- function(count, norm_method = 'tss'){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    ### normalized count matrix
    db.norm <- sweep(count, 1, rowSums(count), FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors
    count_q75 <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    count_N <- rowSums(count)/nrow(count)
    raw_s_factors <- count_q75/count_N
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    ### normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  else if(norm_method == "tmm")
  {
    ## TMM(Trimmed Mean Method)
    ### scale_factors
    count_t <- t(count)
    raw_s_factors <- calcNormFactors(count_t,method = "TMM")
    scale_coeff <- exp((-1/nrow(count)) * sum(log(raw_s_factors)))
    scaled_s_factors <- scale_coeff * raw_s_factors
    # normalized count matrix
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "/")
    count_nor <- db.norm
  }
  
  else if(norm_method == "n-vst") {
    ## Naive Transformation (VST)
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    y <- sqrt( phi * count)
    naive_norm <- log(y + sqrt(1 + y^2))
    total_n_counts <- apply(count, 1, sum)
    log_total_n_counts <- log(total_n_counts)
    db.norm <- apply(naive_norm, 2, function(x){resid(lm(x ~ log_total_n_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(norm_method == "a-vst")
  {
    ## Anscombe's Transformation(VST)
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    ### Took the absolute value of counts as negative under the square root not on Real subspace
    y <- sqrt( abs( (count + (3/8))/((1/phi)-(3/4)) ) )
    anscombe_norm <- log(y + sqrt(1 + y^2))
    total_a_counts <- apply(count, 1, sum)
    log_total_a_counts <- log(total_a_counts)
    db.norm <- apply(anscombe_norm, 2, function(x){resid(lm(x ~ log_total_a_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
  else if(norm_method == "log-vst")
  {
    ## Log Transformation
    varx = apply(count, 2, var)
    meanx = apply(count, 2, mean)
    phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
    log_norm <- log2(count + (1/(2*phi)))
    total_l_counts <- apply(count, 1, sum)
    log_total_l_counts <- log(total_l_counts)
    db.norm <- apply(log_norm, 2, function(x){resid(lm(x ~ log_total_l_counts))} )
    ## All the above normalized counts were negative so reversed their signs
    count_nor <- db.norm
  }
   else if(norm_method == 'none'){
     count_nor <- count
   }
  else
  {
    stop("Please choose a valid normalization method")
  }
  colnames(count_nor) <- colnames(count)
  return(count_nor)
}


cluster.st <- function(count, clustermethod = 'MGC'){
  count_nor <- count
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  # clustering
  dout <- 3
  
  if (clustermethod == 'MGC'){
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      mIQR <- IQR(count_nor[,i])
      list_large <- which(count_nor[,i] > quantile(count_nor[,i], 0.75) + dout*mIQR)
      list_small <- which(count_nor[,i]  == 0)
      
      if (mIQR == 0){list_large <- NULL}
      list_remain <- setdiff(1:sample_num, union(list_large, list_small))
      if (length(list_remain) > 2){
        if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
        {
          count_binary[list_remain,i] <- 1
        }
        else{
          k <- Mclust(count_nor[list_remain,i], G=2, verbose = F)
          if (k$parameters$mean[1]>k$parameters$mean[2])
          {count_binary[list_remain,i] <- 3 - k$classification}
          else {count_binary[list_remain,i] <- k$classification}
        }} else if (length(list_remain) == 1){count_binary[list_remain,i] <- sample(c(1,2), 1)}
      else if (length(list_remain) == 2){ count_binary[list_remain,i] <- order(count_nor[list_remain,i])}
      count_binary[list_large,i] <- 2
      count_binary[list_small,i] <- 1
      
      # If cannot split into two groups, we assume that there is no spatial pattern
      if (sum(count_binary[,i]) == 2*sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
      if (sum(count_binary[,i]) == sample_num)
      {
        count_binary[,i] <- sample(c(1,2),sample_num, replace = TRUE)
        print(i)
      }
    }
  } else if (clustermethod == 'Kmeans'){
    
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      mIQR <- IQR(count_nor[,i])
      list_large <- which(count_nor[,i] > quantile(count_nor[,i], 0.75) + dout*mIQR)
      list_small <- NULL
      
      if (mIQR == 0){list_large <- NULL}
      list_remain <- setdiff(1:sample_num, list_large)
      
      if (length(list_remain)>2){
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else{
        k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
        if (k$centers[1]>k$centers[2])
        {count_binary[list_remain,i] <- 3 - k$cluster}
        else {count_binary[list_remain,i] <- k$cluster}}} else if (length(list_remain) == 1){
        count_binary[list_remain,i] <- sample(c(1,2), 1)
      } else if (length(list_remain) == 1) {count_binary[list_remain,i] <- order(count_nor[list_remain,i])}
      count_binary[list_large,i] <- 2
      
      count_binary[list_small,i] <- 1
    
    if (max(count_binary[,i]) == 1 )
    {
      count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 2
    }
    else if (min(count_binary[,i]) == 2 )
    {
      count_binary[sample(1:sample_num, floor(sample_num/2)), i] <- 1
    }
  }}
  else {stop('Please choose a valid clustering method')}
  colnames(count_binary) <- colnames(count)
  return(count_binary)
}

# run BOOST-Ising 
run.ising <- function(count,sample_info, chain = 1){
  if(is.vector(count) == FALSE){
    stop('Please convert the count matrix for one gene as vector')
  }
  if(length(count) != dim(sample_info)[1]){
    stop('Length of count does not match the number of sample points')
  }
  count_binary <- count
  sample_info <- round(sample_info)
  gene_name <- colnames(count)
  sample_num <- nrow(count)
  # transform to array
  P_nor <- transform_data(count_binary,sample_info)
  # estimate parameters in Ising model
  mean1 = 1
  sigma1 = 2.5
  is_normal = 1
  sigma = 1

  theta_est <- matrix(0, nrow = 5000, ncol = chain)
  omega_est <- matrix(0, nrow = 5000, ncol = chain)
  for (i in 1:chain){
    print(paste0('Chain', i))
    res <- potts_2_omega(P_nor, sigma, mean1, sigma1, rnorm(1, mean = 0, sd = 0.333), rnorm(1, mean = 1, sd = 2.5));
  if (dim(res$theta)[2] > 0 & dim(res$omega)[2] > 0){
      theta_est[,i] <- res$theta[5001:10000,1]
      omega_est[,i] <- res$omega[5001:10000,1]
    }}
  return(list(theta_est, omega_est))
}


# numerically summarize BOOST-Ising result
numeric.sum <- function(ising_result){
theta_est <- ising_result[[1]]
omega_est <- ising_result[[2]]
chain  <- ncol(theta_est)

omega_mean <- mean(omega_est)
theta_mean <- mean(theta_est)

theta_CI_low <- quantile(theta_est,0.025)
theta_CI_high <- quantile(theta_est,0.975)

omega_CI_low <- quantile(theta_est,0.025)
omega_CI_high <- quantile(theta_est,0.975)

theta_sd <- sd(theta_est)
omega_sd <- sd(omega_est)

result <- data.frame(mean = c(theta_mean, omega_mean), CI_low = c(theta_CI_low, omega_CI_low), CI_high = c(theta_CI_high, omega_CI_high), SD = c(theta_sd, omega_sd))
rownames(result) <- c('theta', 'omega')


pvalues_neg <- sum(theta_est>=0)/(5000*chain)
pvalues_pos <- sum(theta_est<=0)/(5000*chain)

BF_neg <- (1 - pvalues_neg)/(pvalues_neg + 0.000000001)
BF_pos <- (1 - pvalues_pos)/(pvalues_pos + 0.000000001)

infer_result <- data.frame(pvalues.neg  = pvalues_neg, BF.neg = BF_neg, pvalues.pos  = pvalues_pos, BF.pos = BF_pos)

return(list(statistics = result, inference = infer_result))
}
  
# graphically summarize BOOST-Ising result
graphical.sum <- function(ising_result){
  theta_est <- ising_result[[1]]
  omega_est <- ising_result[[2]]
  chain  <- ncol(theta_est)
  
  par(mfrow=c(2,2))    # set the plotting area into a 1*2 array
  hist(theta_est, breaks = 50, main = 'Histogram of theta', xlab = 'theta')
  hist(omega_est, breaks = 50, main = 'Histogram of omega',xlab = 'omega')
  
  plot(1:5000, theta_est[, 1], main = 'Traceplot of theta', xlab = 'iteration', col=2, type = 'l', ylab = 'theta')
  if (chain >= 2){
    for (i in 2:chain){
      lines(1:5000, theta_est[,i], col = i +1)
    }
  }
  
  plot(1:5000, omega_est[, 1], main = 'Traceplot of omega', xlab = 'iteration', col=2, type = 'l', ylab = 'omega')
  if (chain >= 2){
    for (i in 2:chain){
      lines(1:5000, omega_est[,i], col = i+1)
    }
  }
}

