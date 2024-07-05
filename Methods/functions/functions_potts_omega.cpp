
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double hamiltonian_c(arma::mat P, NumericMatrix Theta, int H_min, int H_max, int L_min, int L_max);
double hamiltonian2_c_omega(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta, NumericVector omega);
int imax(int a, int b);
int imin(int a, int b);
double dmax(double a, double b);
double dmin(double a, double b);
double myabs(double a);
int getIndex(NumericMatrix potts_info, int h, int l);
int getSum(IntegerMatrix Delta, int k);
IntegerMatrix rectangle2matrix(int H, int L, int H_min, int H_max, int L_min, int L_max, bool rev);
IntegerMatrix potts2ising(IntegerMatrix Delta, int k);
arma::mat potts_c_omega(arma::mat P, NumericMatrix Theta, NumericVector omega);
arma::mat potts2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta);
Rcpp::List potts_2_omega(arma::mat P, double sigma, double omega_mean, double omega_sigma, double theta_initial, double omega_initial);


// [[Rcpp::export]]
Rcpp::List potts_2_omega_laplace(arma::mat P, double sigma,  double omega_mean, double omega_sigma, double theta_initial, double omega_initial) {
  // Read data information
  int H = P.n_rows;
  int L = P.n_cols;
  int Q = P.max();
  
  // Set hyperparameters
  double mu = 0.0;
  //double sigma = 10.0;
  double mu_omega = omega_mean;
  double sigma_omega = omega_sigma;
  
  // Set algorithm settings
  int iter = 10000;
  int burn = iter/2;
  double Theta_s = theta_initial;
  int M = 3;
  
  int i, q, qq, qqq, qqqq, l, h, m, count;
  int count_2 = 10;
  
  double tau = 0.25;
  double hastings = 0;
  double accept = 0;
  double accept_omega = 0;
  IntegerMatrix Delta(H, L);
  NumericMatrix theta_store(iter, Q*(Q - 1)/2);
  NumericMatrix Theta(Q, Q);
  NumericMatrix Theta_temp(Q, Q);
  NumericMatrix Theta_0(Q, Q);
  NumericVector omega(Q);
  NumericVector omega_temp(Q);
  NumericMatrix omega_store(iter, Q);
  arma::mat P_temp(H, L);
  
  
  // Initialization
  for(q = 0; q < Q; q++)
  {
    omega(q) = omega_initial;
    for (qq = 0; qq < Q; qq++)
    {
      Theta(q, qq) = Theta_s;
    }
  }
  for (h = 0; h < H; h++)
  {
    for (l = 0; l < L; l++)
    {
      Delta(h, l) = 1;
    }
  }
  
  omega(Q-1) = 1.0;
  
  // MCMC
  for(i = 0; i < iter; i++)
  {
    // Update omega
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = 0; qq < Q; qq++)
      {
        omega_temp(qq) = omega(qq);
      }
      omega_temp(q) = rnorm(1, omega(q), tau)(0);
      
      for (h = 0; h < H; h++)
      {
        for (l = 0; l < L; l++)
        {
          P_temp(h, l) = P(h, l);
        }
      }
      for (m = 0; m < M; m++)
      {
        P_temp = potts_c_omega(P_temp, Theta, omega_temp);
      }
      hastings = hamiltonian2_c_omega(P_temp, Theta, Delta, omega) - hamiltonian2_c_omega(P, Theta, Delta, omega) + hamiltonian2_c_omega(P, Theta, Delta, omega_temp) - hamiltonian2_c_omega(P_temp, Theta, Delta, omega_temp);
      hastings = hastings - ((omega_temp(q) - mu_omega)*(omega_temp(q) - mu_omega)/2/sigma_omega/sigma_omega - (omega(q) - mu_omega)*(omega(q) - mu_omega)/2/sigma_omega/sigma_omega);
      if (hastings >= log(double(rand()%10001)/10000))
      {
        omega(q) = omega_temp(q);
        if (i > burn) {
          accept_omega++;
        }
      }
    }
    
    // Update Theta
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        for(qqq = 0; qqq < Q; qqq++)
        {
          for (qqqq = 0; qqqq < Q; qqqq++)
          {
            Theta_temp(qqq, qqqq) = Theta(qqq, qqqq);
          }
        }
        Theta_temp(q, qq) = rnorm(1, Theta(q, qq), tau)(0);
        Theta_temp(qq, q) = Theta_temp(q, qq);
        for (h = 0; h < H; h++)
        {
          for (l = 0; l < L; l++)
          {
            P_temp(h, l) = P(h, l);
          }
        }
        for (m = 0; m < M; m++)
        {
          P_temp = potts_c_omega(P_temp, Theta_temp, omega);
        }
        hastings = hamiltonian2_c_omega(P_temp, Theta, Delta, omega) - hamiltonian2_c_omega(P, Theta, Delta,omega) + hamiltonian2_c_omega(P, Theta_temp, Delta,omega) - hamiltonian2_c_omega(P_temp, Theta_temp, Delta, omega);
        hastings = hastings - myabs(Theta_temp(q, qq) - mu)/sigma + myabs(Theta(q, qq) - mu)/sigma;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          Theta(q, qq) = Theta_temp(q, qq);
          Theta(qq, q) = Theta(q, qq);
          if (i > burn) {
            accept++;
          }
        }
      }
    }
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      //Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        theta_store(i, count) = Theta(q, qq);
        count++;
      }
      omega_store(i,q) = omega(q);
    }
    omega_store(i,(Q-1)) = omega(Q-1);
  }
  accept = accept/(iter - burn)/(Q*(Q - 1)/2);
  accept_omega = accept_omega/(iter - burn)/(Q-1);
  return Rcpp::List::create(Rcpp::Named("theta") = theta_store,Rcpp::Named("omega") = omega_store, Rcpp::Named("accept") = accept, Rcpp::Named("accept_omega") = accept_omega);
}

// [[Rcpp::export]]
Rcpp::List potts_2_omega(arma::mat P, double sigma, double omega_mean, double omega_sigma, double theta_initial, double omega_initial) {
  // Read data information
  int H = P.n_rows;
  int L = P.n_cols;
  int Q = P.max();
  
  // Set hyperparameters
  double mu = 0.0;
  // double sigma = 10.0;
  double mu_omega = omega_mean;
  double sigma_omega = omega_sigma;
  
  // Set algorithm settings
  int iter = 10000;
  int burn = iter/2;
  double Theta_s = theta_initial;
  int M = 3;
  
  int i, q, qq, qqq, qqqq, l, h, m, count;
  int count_2 = 10;
  
  double tau = 0.25;
  double hastings = 0;
  double accept = 0;
  double accept_omega = 0;
  IntegerMatrix Delta(H, L);
  NumericMatrix theta_store(iter, Q*(Q - 1)/2);
  NumericMatrix Theta(Q, Q);
  NumericMatrix Theta_temp(Q, Q);
  NumericMatrix Theta_0(Q, Q);
  NumericVector omega(Q);
  NumericVector omega_temp(Q);
  NumericMatrix omega_store(iter, Q);
  arma::mat P_temp(H, L);
  
  
  // Initialization
  for(q = 0; q < Q; q++)
  {
    omega(q) = omega_initial;
    for (qq = 0; qq < Q; qq++)
    {
      Theta(q, qq) = Theta_s;
    }
  }
  for (h = 0; h < H; h++)
  {
    for (l = 0; l < L; l++)
    {
      Delta(h, l) = 1;
    }
  }
  omega(Q-1) = 1.0;
  
  
  // MCMC
  for(i = 0; i < iter; i++)
  {
    // Update omega
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = 0; qq < Q; qq++)
      {
        omega_temp(qq) = omega(qq);
      }
      omega_temp(q) = rnorm(1, omega(q), tau)(0);
      
     for (h = 0; h < H; h++)
     {
       for (l = 0; l < L; l++)
       {
         P_temp(h, l) = P(h, l);
       }
     }
     for (m = 0; m < M; m++)
     {
       P_temp = potts_c_omega(P_temp, Theta, omega_temp);
     }
      hastings = hamiltonian2_c_omega(P_temp, Theta, Delta, omega) - hamiltonian2_c_omega(P, Theta, Delta, omega) + hamiltonian2_c_omega(P, Theta, Delta, omega_temp) - hamiltonian2_c_omega(P_temp, Theta, Delta, omega_temp);
      hastings = hastings - ((omega_temp(q) - mu_omega)*(omega_temp(q) - mu_omega)/2/sigma_omega/sigma_omega - (omega(q) - mu_omega)*(omega(q) - mu_omega)/2/sigma_omega/sigma_omega);
      if (hastings >= log(double(rand()%10001)/10000))
      {
        
        omega(q) = omega_temp(q);
        if (i > burn) {
          accept_omega++;
        }
      }
    }
    
    // Update Theta
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        for(qqq = 0; qqq < Q; qqq++)
        {
          for (qqqq = 0; qqqq < Q; qqqq++)
          {
            Theta_temp(qqq, qqqq) = Theta(qqq, qqqq);
          }
        }
        Theta_temp(q, qq) = rnorm(1, Theta(q, qq), tau)(0);
        Theta_temp(qq, q) = Theta_temp(q, qq);
        for (h = 0; h < H; h++)
        {
          for (l = 0; l < L; l++)
          {
            P_temp(h, l) = P(h, l);
          }
        }
        for (m = 0; m < M; m++)
        {
          P_temp = potts_c_omega(P_temp, Theta_temp, omega);
        }
        hastings = hamiltonian2_c_omega(P_temp, Theta, Delta, omega) - hamiltonian2_c_omega(P, Theta, Delta,omega) + hamiltonian2_c_omega(P, Theta_temp, Delta,omega) - hamiltonian2_c_omega(P_temp, Theta_temp, Delta, omega);
        hastings = hastings - (Theta_temp(q, qq) - mu)*(Theta_temp(q, qq) - mu)/2/sigma/sigma + (Theta(q, qq) - mu)*(Theta(q, qq) - mu)/2/sigma/sigma;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          Theta(q, qq) = Theta_temp(q, qq);
          Theta(qq, q) = Theta(q, qq);
          if (i > burn) {
            accept++;
          }
        }
      }
    }
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      //Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < Q - 1; q++)
    {
      for(qq = q + 1; qq < Q; qq++)
      {
        theta_store(i, count) = Theta(q, qq);
        count++;
      }
      omega_store(i,q) = omega(q);
    }
    omega_store(i,(Q-1)) = omega(Q-1);
  }
  accept = accept/(iter - burn)/(Q*(Q - 1)/2);
  accept_omega = accept_omega/(iter - burn)/(Q-1);
  return Rcpp::List::create(Rcpp::Named("theta") = theta_store,Rcpp::Named("omega") = omega_store, Rcpp::Named("accept") = accept, Rcpp::Named("accept_omega") = accept_omega);
}

// [[Rcpp::export]]
IntegerMatrix rectangle2matrix(int H, int L, int H_min, int H_max, int L_min, int L_max, bool rev) {
  IntegerMatrix Delta(H, L);
  int h, l;
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if (h >= H_min - 1 & h < H_max & l >= L_min - 1 & l < L_max)
      {
        if(rev) 
        {
          Delta(h, l) = 0;
        }
        else
        {
          Delta(h, l) = 1;
        }
      }
      else
      {
        if(rev) 
        {
          Delta(h, l) = 1;
        }
        else
        {
          Delta(h, l) = 0;
        }
      }
    }
  }
  return Delta;
}

// [[Rcpp::export]]
double hamiltonian_c(arma::mat P, NumericMatrix Theta, int H_min, int H_max, int L_min, int L_max) {
  double hamiltonian = 0;
  //int L = P.n_cols;
  //int H = P.n_rows;
  int h, l;
  for(h = H_min - 1; h < H_max - 1; h++)
  {  
    for(l = L_min - 1; l < L_max - 1; l++)
    {
      if(P(h, l) != 0)
      {
        if(P(h, l) != P(h, l + 1))
        {
          hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
        }
        if(P(h, l) != P(h + 1, l))
        {
          hamiltonian = hamiltonian + Theta(P(h + 1, l) - 1, P(h, l) - 1);
        }
      }
    }
  }
  h = H_max - 1;
  for(l = L_min - 1; l < L_max - 1; l++)
  {
    if(P(h, l) != 0 && P(h, l) != P(h, l + 1))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
    }
  }
  l = L_max - 1;
  for(h = H_min - 1; h < H_max - 1; h++)
  {
    if(P(h, l) != 0 && P(h, l) != P(h + 1, l))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h + 1, l) - 1);
    }
  }
  return (hamiltonian);
}

// [[Rcpp::export]]
double hamiltonian2_c_omega(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta, NumericVector omega) {
  double hamiltonian = 0;
  int L = P.n_cols;
  int H = P.n_rows;
  int h, l;
  for(h = 0; h < H - 1; h++)
  {  
    for(l = 0; l < L - 1; l++)
    {
      if(P(h, l) != 0 && Delta(h, l) != 0)
      {
        if(Delta(h, l + 1) != 0 && P(h, l) != P(h, l + 1))
        {
          hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
        }
        if(Delta(h + 1, l) != 0 && P(h, l) != P(h + 1, l))
        {
          hamiltonian = hamiltonian + Theta(P(h + 1, l) - 1, P(h, l) - 1);
        }
      }
    }
  }
  h = H - 1;
  for(l = 0; l < L - 1; l++)
  {
    
    if(P(h, l) != 0 && Delta(h, l) != 0 && Delta(h, l + 1) != 0 && P(h, l) != P(h, l + 1))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h, l + 1) - 1);
    }
  }
  l = L - 1;
  for(h = 0; h < H - 1; h++)
  {
    if(P(h, l) != 0 && Delta(h, l) != 0 && Delta(h + 1, l) !=0 && P(h, l) != P(h + 1, l))
    {
      hamiltonian = hamiltonian + Theta(P(h, l) - 1, P(h + 1, l) - 1);
    }
  }
  for(h = 0; h < H; h++)
  {  
    for(l = 0; l < L; l++)
    {
      if(P(h, l) != 0 && Delta(h, l) != 0)
      {
        hamiltonian = hamiltonian + omega(P(h,l) - 1);
      }
    }
  }
  return (hamiltonian);
}

// [[Rcpp::export]]
arma::mat potts_c_omega(arma::mat P, NumericMatrix Theta, NumericVector Omega) {
  int L = P.n_cols;
  int H = P.n_rows;
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int h, l, q;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(P(h, l) != 0)
      {
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = Omega(q);
          if(l < L - 1 && (P(h, l + 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l + 1) - 1, q);
          }
          if(l > 0 && (P(h, l - 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l - 1) - 1, q);
          }
          if(h < H - 1 && (P(h + 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h + 1, l) - 1, q);
          }
          if(h > 0 && (P(h - 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h - 1, l) - 1, q);
          }
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat potts2_c(arma::mat P, NumericMatrix Theta, IntegerMatrix Delta) {
  int L = P.n_cols;
  int H = P.n_rows;
  int Q = Theta.nrow();
  NumericVector prob_temp(Q);
  IntegerVector state(Q);
  int h, l, q;
  double temp = 0;
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(P(h, l) != 0 && Delta(h, l) != 0)
      {
        for(q = 0; q < Q; q++)
        {
          prob_temp(q) = 0;
          if(l < L - 1 && Delta(h, l + 1) != 0 && (P(h, l + 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l + 1) - 1, q);
          }
          if(l > 0 && Delta(h, l - 1) != 0 && (P(h, l - 1) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h, l - 1) - 1, q);
          }
          if(h < H - 1 && Delta(h + 1, l) != 0 && (P(h + 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h + 1, l) - 1, q);
          }
          if(h > 0 && Delta(h - 1, l) != 0 && (P(h - 1, l) - 1) != q)
          {
            prob_temp(q) = prob_temp(q) + Theta(P(h - 1, l) - 1, q);
          }
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat MRFi(arma::mat P, int H, int L, double e, double f, int iter) {
  int h, l, i;
  double prob;
  for (i = 0; i < iter; i++)
  {
    for (h = 0; h < H; h++)
    {
      for (l = 0; l < L; l++)
      {
        P(h, l) = 0;
        prob = exp(e + f*(P(imin(h + 1, H - 1), l) + P(imax(h - 1, 0), l) + P(h, imin(l + 1, L - 1)) + P(h, imax(l - 1, 0))));
        prob = prob/(1 + prob);
        P(h, l) = rbinom(1, 1, prob)(0);
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
arma::mat MRFp(arma::mat P, int H, int L, int Q, double f, int iter) {
  int h, l, q, i;
  double temp;
  IntegerVector state(Q);
  NumericVector prob_temp(Q);
  for(q = 0; q < Q; q++)
  {
    state(q) = q;
  }
  for (i = 0; i < iter; i++)
  {
    for (h = 0; h < H; h++)
    {
      for (l = 0; l < L; l++)
      {
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = 0;
        }
        if (l < L - 1) {
          prob_temp(P(h, l + 1) - 1) = prob_temp(P(h, l + 1) - 1) + f;
        }
        if (l > 0) {
          prob_temp(P(h, l - 1) - 1) = prob_temp(P(h, l - 1) - 1) + f;
        }
        if (h < H - 1) {
          prob_temp(P(h + 1, l) - 1) = prob_temp(P(h + 1, l) - 1) + f;
        }
        if (h > 0) {
          prob_temp(P(h - 1, l) - 1) = prob_temp(P(h - 1, l) - 1) + f;
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = exp(prob_temp(q));
        }
        temp = 0;
        for (q = 0; q < Q; q++)
        {
          temp = temp + prob_temp(q);
        }
        for (q = 0; q < Q; q++)
        {
          prob_temp(q) = prob_temp(q)/temp;
        }
        P(h, l) = RcppArmadillo::sample(state, 1, true, prob_temp)(0) + 1;
      }
    }
  }
  return P;
}

// [[Rcpp::export]]
IntegerMatrix potts2ising(IntegerMatrix Delta, int k) {
  int L = Delta.ncol();
  int H = Delta.nrow();
  int l, h;
  IntegerMatrix Delta_k(H, L);
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(Delta(h, l) == k)
      {
        Delta_k(h, l) = 1;
      }
      else
      {
        Delta_k(h, l) = 0;
      }
    }
  }
  return Delta_k;
}

// [[Rcpp::export]]
int getSum(IntegerMatrix Delta, int k) {
  int L = Delta.ncol();
  int H = Delta.nrow();
  int sum = 0;
  int h, l;
  for(h = 0; h < H; h++)
  {
    for(l = 0; l < L; l++)
    {
      if(Delta(h, l) == k)
      {
        sum++;
      }
    }
  }
  return sum;
}

// [[Rcpp::export]]
int getIndex(NumericMatrix potts_info, int h, int l) {
  int index = 0;
  while(h != potts_info(index, 0) || l != potts_info(index, 1))
  {
    index = index + 1;
  }
  return index + 1;
}

// [[Rcpp::export]]
int imax(int a, int b) {
  if(a > b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
int imin(int a, int b) {
  if(a < b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
double dmax(double a, double b) {
  if(a > b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
double dmin(double a, double b) {
  if(a < b) 
  {
    return a;
  }
  else
  {
    return b;
  }
}

// [[Rcpp::export]]
double myabs(double a) {
  if(a < 0) 
  {
    return -a;
  }
  else
  {
    return a;
  }
}
