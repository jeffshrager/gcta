data {

  // Parameter counts
  int<lower=0> N_bm;  // number of biomarkers (with treatment interaction effects)
  int<lower=0> N_tx;  // number of treatments types
  
  // hyperparameters
  // random effects hyperparameters
  real<lower=0> eta_u;            // correlation hyperparameter
  real loc_mu_u1;                 // location parameters for slope mean
  real scale_mu_u1;               // scale parameters for slope mean
  vector<lower=0>[2] scale_tau_u; // scale parameter for intercept/slope std

  // fixed effects hyperparameters
  vector[N_bm] loc_bm;            // location of biomarker effects
  matrix[N_bm, N_bm] scale_bm;    // scale of biomarker effects
  vector[N_tx] loc_tx;            // location of treatment type effects
  matrix[N_tx, N_tx] scale_tx;    // scale of treatment type effects
  vector[N_bm * N_tx] loc_int;    // location of interaction effects
  matrix[N_bm * N_tx,
	 N_bm * N_tx] scale_int;  // scale of interaction effects

  // measurement uncertainty hyperparameter
  real scale_sigma_eps;           // scale of measurement error term
}


transformed data {
  // number of treatment-biomarker interactions
  int<lower=0> N_int = N_bm * N_tx;    
}


parameters {
  
  // fixed (biomarker/tx-level) effects
  real<lower=0, upper=100> beta0;  // intercept
  vector[N_bm] beta_bm;            // biomarker effects
  vector[N_tx] beta_tx;            // treatment effects
  vector[N_bm * N_tx] beta_int;    // interaction effects
  
  // random (patient-level) effects hyperparameters
  real mu_u1;                  // mean of slope
  real<lower=0, upper=1> b_u;  // rho_u / 2 + 1
  vector<lower=0>[2] tau_u;    // scale parameters for intercept and slope
  
  // measurement standard deviation  
  real<lower=0> sigma_eps;
  
}


transformed parameters {
  // random effects
  real<lower=-1, upper=1> rho_u;
  rho_u = 2 * b_u - 1;
}

model {
  
  //----------------------------------------------------------------------------
  // priors
  //----------------------------------------------------------------------------

  // random effects mean of intercept and slope
  mu_u1 ~ normal(loc_mu_u1, scale_mu_u1);

  // covariance between random effects intercept and slope
  tau_u ~ cauchy(0, scale_tau_u);
  b_u ~ beta(eta_u, eta_u);

  // fixed effects predictor coefficients
  // beta0 is uniform between 0 and 100
  beta_bm ~ multi_normal(loc_bm, scale_bm);
  beta_tx ~ multi_normal(loc_tx, scale_tx);
  beta_int ~ multi_normal(loc_int, scale_int);

  // measurement error
  sigma_eps ~ cauchy(0, scale_sigma_eps);
}
