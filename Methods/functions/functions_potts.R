
array2matrix_r = function(theta) {
  Q <- floor(sqrt(2*length(theta))) + 1;
  Theta <- matrix(0, nrow = Q, ncol = Q);
  count <- 1;
  for (q in 1:(Q - 1)) {
    for (qq in (q + 1):Q) {
      Theta[q, qq] <- theta[count];
      Theta[qq, q] <- theta[count];
      count <- count + 1;
    }
  }
  return (Theta);
}

bayes_fdr <- function(PPI, alpha) {
  PPI_sorted <- sort(PPI, decreasing = TRUE);
  k <- 1;
  fdr <- 0;
  while(fdr < alpha) {
    fdr <- mean(1 - PPI_sorted[1:k]);
    k <- k + 1;
  }
  return(PPI_sorted[k]);
}