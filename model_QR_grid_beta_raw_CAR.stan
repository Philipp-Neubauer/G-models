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
  real sparse_car_lpdf(vector phi, real tau, real alpha, 
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
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phi - alpha * (phit_W * phi)));
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
  int<lower=0> W_n;
  matrix<lower=0, upper=1>[N_grid,N_grid] W;
 }
transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  vector[N_grid] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N_grid] D_sparse;     // diagonal of D (number of neigbors for each site)

  // thin and scale the QR decomposition
  Q_ast = qr_Q(COVS)[, 1:K] * sqrt(N - 1.0);
  R_ast = qr_R(COVS)[1:K, ] / sqrt(N - 1.0);
  R_ast_inverse = inverse(R_ast);
  
   { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N_grid - 1)) {
      for (j in (i + 1):N_grid) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N_grid) D_sparse[i] = sum(W[i]);
  {
    vector[N_grid] invsqrtD;  
    for (i in 1:N_grid) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form_diag(W,invsqrtD));
  }
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
  real<lower=0> grid_sig;
  vector[N_grid] grid_fx;
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  real omean;
  real zmean;
  vector[K] theta;
}
transformed parameters{
 vector[N] imean;
 vector<lower=0>[N_notzero] mu;
 vector<lower=0>[N_notzero] a;
 vector<lower=0>[N_notzero] b;
  
 imean = Q_ast*theta+ island_sig*island_fx[ISLAND]+group_sig*group_fx[GROUP]+ grid_fx[GRID];

 mu = inv_logit(omean+imean[ii_notzero]);
 a  =    mu  .* iscale[ISLAND[ii_notzero]];
 b  = (1-mu) .* iscale[ISLAND[ii_notzero]];
 

}
model {
  
    // model for within grid sd, hierarchical across islands
  iscale ~ cauchy(0,gscale);
  gscale ~ cauchy(0,5);
  
  IMEAN ~ beta(a,b);
  
  //hurdle model for single measurement outcomes
  ZEROONE ~ bernoulli_logit(zmean+scale*imean[ii_zeroone]);
  
  //regression coeffs
  omean ~ cauchy(0,10);
  zmean ~ cauchy(0,10);
  theta ~ cauchy(0,5);
  scale ~ cauchy(1,2);
  
  island_fx ~ normal(0,1);
  island_sig ~ cauchy(0,1);
  group_fx ~ normal(0,1);
  group_sig ~ cauchy(0,1);
  grid_sig ~ cauchy(0,1);
  grid_fx ~ sparse_car(grid_sig, rho, W_sparse, D_sparse, lambda, N_grid, W_n);
  rho ~ beta(1, 1);
  
}
generated quantities{
  //vector[N] log_lik;
  vector[K] beta;
  vector[N_notzero] resid_y;
  beta = R_ast_inverse * theta; // coefficients
  resid_y = logit(IMEAN) - (omean+imean[ii_notzero]);
  // beta_pdf for means and draws
  //for (i in 1:N) log_lik[i] =+ bernoulli_logit_lpmf(ZEROONE[i] | scale*imean[i]);
  //for (i in 1:N_notzero) log_lik[ii_notzero[i]] = beta_lpdf(IMEAN[i] | a[i], b[i]);
  // Hurdle part of the likelihood for draws
}
