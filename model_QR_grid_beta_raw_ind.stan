data {
  int<lower = 0> N_notzero;
  int<lower = 0> N_groups;
  int<lower = 0> N_islands;
  int<lower = 0> N_grid;
  int<lower = 0> N;
  int<lower = 0> K;
  vector<lower = 0, upper = 1>[N_notzero] IMEAN;
  int<lower = 0, upper = 1> ZEROONE[N];
  int<lower = 0> ii_notzero[N_notzero];
  matrix[N, K] COVS;
  int<lower = 1, upper = N_islands> ISLAND[N];
  int<lower = 1, upper = N_groups> GROUP[N];
  int<lower = 1, upper = N> GRID[N];
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
  real zmean;
  real scale;
  vector[N_groups] group_fx;
  real<lower=0> island_sig;
  real<lower=0> group_sig;
  vector[N_grid] grid_fx;
  vector<lower=0>[N_islands] grid_sigs;
  real omean;
  matrix[N_islands,K] theta;
}
transformed parameters{
 vector[N] imean;
 vector<lower=0>[N_notzero] mu;
 vector<lower=0>[N_notzero] a;
 vector<lower=0>[N_notzero] b;
  
 imean = rows_dot_product(Q_ast,theta[ISLAND,1:K])+island_sig*island_fx[ISLAND]+group_sig*group_fx[GROUP]+grid_sigs[ISLAND].*grid_fx[GRID];

 mu = inv_logit(omean+imean[ii_notzero]);
 a  =    mu  .* iscale[ISLAND[ii_notzero]];
 b  = (1-mu) .* iscale[ISLAND[ii_notzero]];
}
model {
  
    // model for within grid sd, hierarchical across islands
  
  iscale ~ cauchy(0,10);
  
  IMEAN ~ beta(a,b);
  
  //hurdle model for single measurement outcomes
  ZEROONE ~ bernoulli_logit(zmean+scale*imean);
  
  //regression coeffs
  omean ~ cauchy(0,5);
  zmean ~ cauchy(0,5);
  scale ~ cauchy(0,5);
  
  for (i in 1:N_islands) theta[i,1:K] ~ cauchy(0,10);
  
  island_fx ~ normal(0,1);
  island_sig ~ cauchy(0,1);
  group_fx ~ normal(0,1);
  group_sig ~ cauchy(0,1);
  grid_fx ~ normal(0,1);
  grid_sigs ~ cauchy(0,1);
  
}
generated quantities{
  vector[N] log_lik;
  //vector[N] resid_z;
  vector[N_notzero] resid_y;
  vector[K] beta[N_islands];
    for (i in 1:N_islands) beta[i] = R_ast_inverse*to_vector(theta[i,1:K]); // coefficients
  // bernoulli for 0-1 
  for (i in 1:N) log_lik[i] =+ bernoulli_logit_lpmf(ZEROONE[i] | zmean+scale*imean[i]);
  // beta_pdf for positive draws
  for (i in 1:N_notzero) log_lik[ii_notzero[i]] = beta_lpdf(IMEAN[i] | a[i], b[i]);
  resid_y = logit(IMEAN) - (omean+imean[ii_notzero]);
   // Hurdle part of the likelihood for draws
}
