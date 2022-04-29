data {

  // Parameter counts
  int<lower=0> N_bm;  // number of biomarkers (with treatment interaction effects)
  int<lower=0> N_tx;  // number of treatments types
    
  // hyperparameters
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
  real beta0;                      // intercept
  vector[N_bm] beta_bm;            // biomarker effects
  vector[N_tx] beta_tx;            // treatment effects
  vector[N_bm * N_tx] beta_int;    // interaction effects
  
  // measurement standard deviation  
  real<lower=0> sigma_eps;
  
}


model {

  //----------------------------------------------------------------------------
  // priors
  //----------------------------------------------------------------------------

  // fixed effects predictor coefficients
  beta0 ~ normal(0, 10);
  beta_bm ~ multi_normal(loc_bm, scale_bm);
  beta_tx ~ multi_normal(loc_tx, scale_tx);
  beta_int ~ multi_normal(loc_int, scale_int);

  // measurement error
  sigma_eps ~ cauchy(0, scale_sigma_eps);
}
