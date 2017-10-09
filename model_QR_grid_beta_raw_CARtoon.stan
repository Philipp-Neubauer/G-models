functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real alpha, 
    int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (sum(ldet_terms)
                    - (phit_D * phi - alpha * (phit_W * phi)));
  }
}
data {
  int<lower = 0> N_notzero;
  int<lower = 0> N_zeroone;
  int<lower = 0> N_groups;
  int<lower = 0> N_islands;
  int<lower = 0> N_grid;
  int<lower = 0> N;
  int<lower = 0> K;
  vector<lower = 0, upper = 1>[N_notzero] IMEAN;
  int<lower = 0, upper = 1> ZEROONE[N_zeroone];
  int<lower = 0> ii_notzero[N_notzero];
  int<lower = 0> ii_map[N_notzero];
  int<lower = 0> ii_zeroone[N_zeroone];
  matrix[N, K] COVS;
  int<lower = 1, upper = N_islands> ISLAND[N];
  int<lower = 1, upper = N_groups> GROUP[N];
  int<lower = 1, upper = N> GRID[N];
  int<lower=0> W_n[N_islands];
  int<lower=0> W_ns;
  int<lower=0> grids[N_islands];
  vector[N_grid] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  int W_sparse[W_ns, 2];   // adjacency pairs
  vector[N_grid] D_sparse;     // diagonal of D (number of neigbors for each site)
 }
transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  
  // thin and scale the QR decomposition
  Q_ast = qr_Q(COVS)[, 1:K] * sqrt(N - 1.0);
  R_ast = qr_R(COVS)[1:K, ] / sqrt(N - 1.0);
  R_ast_inverse = inverse(R_ast);

}
parameters {
  vector[N_islands] island_fx;
  vector<lower=0>[N_islands] iscale;
  real<lower=0> gscale;
  //real<lower=0> vscale;
  real<lower=0> scale;
  vector[N_groups] group_fx;
  real<lower=0> island_sig;
  real<lower=0> group_sig;
  vector<lower=0>[N_islands] grid_sig;
  real<lower=0> grid_var;
  real<lower=0> grid_mean;
  vector[N_grid] grid_fx;
  real<lower=0, upper=1> rho[N_islands]; 
  real<lower=0, upper=1> rho_mean; 
  real<lower=0> rho_prec;  
  real omean;
  real zmean;
  vector[K] theta;
}
transformed parameters{
 vector[N] imean;
 vector<lower=0>[N_notzero] mu;
 vector<lower=0>[N_notzero] a;
 vector<lower=0>[N_notzero] b;
  
 imean = Q_ast*theta+ island_sig*island_fx[ISLAND]+group_sig*group_fx[GROUP]+ grid_sig[ISLAND].*grid_fx[GRID];

 mu = inv_logit(omean+imean[ii_notzero]);
 a  =    mu  .* iscale[ISLAND[ii_notzero]];
 b  = (1-mu) .* iscale[ISLAND[ii_notzero]];
 

}
model {
  int pos;
  int ppos;
  real a_r;
  real b_r;
  
 
    // model for within grid sd, hierarchical across islands
  iscale ~ cauchy(0,gscale);
  gscale ~ cauchy(0,10);
  
  IMEAN ~ beta(a,b);
  
  //hurdle model for single measurement outcomes
  ZEROONE ~ bernoulli_logit(zmean+scale*imean[ii_zeroone]);
  
  //regression coeffs
  omean ~ cauchy(0,10);
  zmean ~ cauchy(0,10);
  theta ~ cauchy(0,5);
  scale ~ cauchy(1,5);
  
  island_fx ~ normal(0,1);
  island_sig ~ cauchy(0,2);
  group_fx ~ normal(0,1);
  group_sig ~ cauchy(0,2);
  grid_sig ~ cauchy(grid_mean,grid_var);
  grid_mean ~ cauchy(0,2);
  grid_var ~ cauchy(0,2);
  
  pos=1;
  ppos=1;
  for (i in 1:N_islands) {
    segment(grid_fx,pos,grids[i]) ~ sparse_car(rho[i], segment(W_sparse,ppos,W_n[i]), segment(D_sparse,pos,grids[i]), segment(lambda,pos,grids[i]), grids[i], W_n[i]);
    pos = pos + grids[i];
    ppos = ppos + W_n[i];
  }
  
 rho_mean  ~ beta(1, 1);
 rho_prec  ~ cauchy(0,1);
 a_r  = rho_mean  * rho_prec;
 b_r  = (1-rho_mean)  * rho_prec;
  
 rho ~ beta(a_r, b_r);
  
}
generated quantities{
  vector[N_zeroone] log_lik;
  vector[K] beta;
  vector[N_notzero] resid_y;
  beta = R_ast_inverse * theta; // coefficients
  resid_y = logit(IMEAN) - (omean+imean[ii_notzero]);
  
  // Hurdle part of the likelihood for 0-1 draws
  for (i in 1:N_zeroone) log_lik[i] =+ bernoulli_logit_lpmf(ZEROONE[i] | zmean + scale*imean[ii_zeroone[i]]);
   // beta_pdf for means and draws
  for (i in 1:N_notzero) log_lik[ii_map[i]] = beta_lpdf(IMEAN[i] | a[i], b[i]);
}
