# load packages

library(mclust)
library(lattice)
library(edgeR)

source("functions/functions_potts.R");
Rcpp::sourceCpp('functions/functions_potts_omega.cpp');

transform_data <- function(count, loc){
  array_y <- array(0, dim = c((max(loc[,1])-min(loc[,1]))+1,(max(loc[,2])-min(loc[,2]))+1,ncol(count)))
  
   loc1 <- min(loc[,1])
   loc2 <- min(loc[,2])
  #loc1 = 0
  #loc2 = 0
  for (i in 1:ncol(count)){
    for (j in 1:nrow(loc)){
      array_y[(loc[j,1]-loc1)+1,(loc[j,2]-loc2)+1, i] <- count[j,i]
    }
  }
  return(array_y)
}


# filter data
filter_count <- function(count, sample_info, min_total = 10,min_percentage = 0.1){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  if (sum(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num) > 0){
  if(sum(rowSums(count) < min_total) == 0){
    sample_f <- sample_info
    count_f <- count[,-(which(colSums(count == 0) > (1-min_percentage)*sample_num))]
  }
  else{
  sample_f <- sample_info[-which(rowSums(count)<min_total),]
  count_f <- count[-which(rowSums(count)<min_total),-(which(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num))]
  }}
  else{
    if(sum(rowSums(count) < min_total) == 0){
      sample_f <- sample_info
      count_f <- count
    }
    else{
      sample_f <- sample_info[-which(rowSums(count)<min_total),]
      count_f <- count[-which(rowSums(count)<min_total),]
    }
  }
  return(list(sample_f,count_f))
}


# main function
Boost_Ising <- function(count,sample_info, norm_method = 'tss', clustermethod = 'MGC')
{
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
    db.norm <- sweep(count, 1, raw_s_factors, FUN = "*")
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
  else if(norm_method == "log")
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
  else
  {
    stop("Please choose a normalization method")
  }
  
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
      else
      {
      k <- Mclust(count_nor[list_remain,i], G=2, verbose = FALSE)
      if (k$parameters$mean[1]>k$parameters$mean[2])
      {count_binary[list_remain,i] <- 3 - k$classification}
      else {count_binary[list_remain,i] <- k$classification}
    }}
      else if (length(list_remain) == 1){
        count_binary[list_remain,i] <- sample(c(1, 2), 1)
      }
     else if (length(list_remain) == 2){
       count_binary[list_remain,i] <- order(count_nor[list_remain,i])
     }
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
  }
    
  # Kmeans clustering
  
  if (clustermethod == 'Kmeans'){
   
    count_binary <- matrix(0, nrow = sample_num, ncol = gene_num) 
    
    set.seed(123)
    for (i in 1:gene_num){
      mIQR <- IQR(count_nor[,i])
    list_large <- which(count_nor[,i] > quantile(count_nor[,i], 0.75) + dout*mIQR)
    list_small <- NULL

    if (mIQR == 0){list_large <- NULL}
    list_remain <- setdiff(1:sample_num, list_large)
    
   if (length(list_remain) > 2){
      if (min(count_nor[list_remain,i]) == max(count_nor[list_remain,i]))
      {
        count_binary[list_remain,i] <- 1
      }
      else{
      k <- kmeans(count_nor[list_remain,i], 2, nstart = 25)
      if (k$centers[1]>k$centers[2])
      {count_binary[list_remain,i] <- 3 - k$cluster}
      else {count_binary[list_remain,i] <- k$cluster}}}
      else if (length(list_remain) == 1){
        count_binary[list_remain,i] <- sample(c(1, 2), 1)
      }
     else if (length(list_remain) == 2){
       count_binary[list_remain,i] <- order(count_nor[list_remain,i])
     }
      
      count_binary[list_large,i] <- 2
      count_binary[list_small,i] <- 1
    }

    if (max(count_binary[,i]) == 1 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 2
    }
    else if (min(count_binary[,i]) == 2 )
    {
      count_binary[sample(1:sample_num, 100), i] <- 1
    }
  }
  
  # transform to array
  P_nor <- transform_data(count_binary,sample_info)
  
  # estimate parameters in Potts model
  theta_est <- array(0,dim=c(5000, gene_num))
  omega_est <- array(0,dim=c(5000, gene_num))
  
  theta_CI_low <- numeric(gene_num)
  theta_CI_high <- numeric(gene_num)
  
  theta_mean <- numeric(gene_num)
  omega_mean <- numeric(gene_num)
  
  omega_CI_low <- numeric(gene_num)
  omega_CI_high <- numeric(gene_num)
      mean1 = 1
      sigma1 = 2.5
      is_normal = 1
      sigma = 1
  count_2 <- 10
  for (i in 1:gene_num){
   
    res <- potts_2_omega(P_nor[,,i], sigma, mean1, sigma1, rnorm(1,0,0.333), rnorm(1,1,2.5));
 
    if (dim(res$theta)[2] > 0 & dim(res$omega)[2] > 0){
    theta_est[,i] <- res$theta[5001:10000,1]
    omega_est[,i] <- res$omega[5001:10000,1]
    
    theta_CI_low[i] <- quantile(theta_est[,i],0.025)
    theta_CI_high[i] <- quantile(theta_est[,i],0.975)
    
    omega_CI_low[i] <- quantile(omega_est[,i],0.025)
    omega_CI_high[i] <- quantile(omega_est[,i],0.975)

    theta_mean[i] <- mean(res$theta[5001:10000,1])
    omega_mean[i] <- mean(res$omega[5001:10000,1])
    }
   if (floor(i*100/gene_num) == count_2)
  {
    print(paste0(count_2, '% has been done'))
    count_2 = count_2 + 10;
  }
  }
  
  
  # detect SV genes
  pvalues_neg <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_neg[i] <- sum(theta_est[,i]>=0)/5000
  }

  pvalues_pos <- numeric(gene_num)
  for (i in 1:gene_num){
    pvalues_pos[i] <- sum(theta_est[,i]<=0)/5000
  }
  
  BF_neg <- pvalues_pos/(pvalues_neg + 0.000000001)
  BF_pos <- pvalues_neg/(pvalues_pos + 0.000000001)
  results <- data.frame(theta_mean, theta_CI_low, theta_CI_high, omega_mean,omega_CI_low,omega_CI_high, BF_neg, BF_pos)
  rownames(results) <- colnames(count)
  return(results)
}  
