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
  int<lower=0> N_edges;
  real<lower=0> scaling_factor;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]
}
transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  
  // thin and scale the QR decomposition
  Q_ast = qr_Q(COVS)[, 1:K] * sqrt(N - 1);
  R_ast = qr_R(COVS)[1:K, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
  
}
parameters {
  //vector<lower=0, upper=0.5>[N_notzero] esig;
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
  vector[N_grid-1] sp_raw;
  real<lower=0> sigma;        // overall standard deviation
  //real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance
  real omean;
  real zmean;
  vector[K] theta;
}
transformed parameters{
 vector[N] imean;
 //vector[N_grid] convolved_re;
 vector<lower=0>[N_notzero] mu;
 vector<lower=0>[N_notzero] a;
 vector<lower=0>[N_notzero] b;
 vector[N_grid] sp_fx;
 
 sp_fx[1:(N_grid - 1)] = sp_raw;
 sp_fx[N_grid] = -sum(sp_raw);
  
  
 imean = Q_ast*theta+island_sig*island_fx[ISLAND]+group_sig*group_fx[GROUP]+ sp_fx[GRID] * sigma/sqrt(scaling_factor)+ grid_fx[GRID]*grid_sig;

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
  
  target += -0.5 * dot_self(sp_fx[node1] - sp_fx[node2]);
  
  //regression coeffs
  omean ~ cauchy(0,10);
  zmean ~ cauchy(0,10);
  theta ~ cauchy(0,5);
  scale ~ cauchy(1,2);
  
  island_fx ~ normal(0,1);
  island_sig ~ cauchy(0,5);
  group_fx ~ normal(0,1);
  group_sig ~ cauchy(0,5);
  
  grid_fx ~ normal(0,1);
  grid_sig ~ cauchy(0,5);
  sigma ~ cauchy(0,5);
  //rho ~ beta(2, 2);
  
}
generated quantities{
  //vector[N] log_lik;
  vector[K] beta;
  beta = R_ast_inverse * theta; // coefficients
  // beta_pdf for means and draws
  //for (i in 1:N) log_lik[i] =+ bernoulli_logit_lpmf(ZEROONE[i] | scale*imean[i]);
  //for (i in 1:N_notzero) log_lik[ii_notzero[i]] = beta_lpdf(IMEAN[i] | a[i], b[i]);
  // Hurdle part of the likelihood for draws
}
