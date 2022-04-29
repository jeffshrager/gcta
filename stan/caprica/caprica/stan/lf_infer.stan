functions {

  vector mean_res(int N_res,       // number of responses
		  int N_bm,        // number of biomarkers
		  int N_tx,        // number of treatments
		  int N_int,       // number of interactions
		  vector t_res,    // time of response
		  int[] pt_index,  // index of patient for response
		  matrix bm_ind,   // N_pt by N_bm indicator matrix
		  matrix tx_ind,   // N_pt by N_tx indicator matrix
		  matrix int_ind,  // N_pt by N_int indicator matrix
		  real beta0,      // constant fixed effect term
		  vector beta_bm,  // biomarker effects term
		  vector beta_tx,  // treatment effects term
		  vector beta_int) // interaction effects term
  {
    vector[N_res] intercept;
    vector[N_res] fixed_effects_slope;
    
    // temporary variables for modeling each response datapoint
    real t;                 // time of response
    int pt;                 // patient index of response
    
    vector[N_bm] x_bm;      // biomarker predictors
    vector[N_tx] x_tx;      // treatment predictors
    vector[N_int] x_int;    // interaction predictors
    
    // build N_res-sized intercept and slope terms
    for (i in 1:N_res) {
    
      t = t_res[i];
      pt = pt_index[i];
    
      intercept[i] = beta0;

      // build predictors
      x_bm = (bm_ind[pt, :] * t)';
      x_tx = (tx_ind[pt, :] * t)';
      x_int = (int_ind[pt, :] * t)';
    
      fixed_effects_slope[i] =
	+ dot_product(beta_bm, x_bm)
	+ dot_product(beta_tx, x_tx)
	+ dot_product(beta_int, x_int);
    }
    return intercept + fixed_effects_slope;
  }
}


data {

  // Parameter counts
  int<lower=0> N_bm;  // number of biomarkers (with treatment interaction effects)
  int<lower=0> N_tx;  // number of treatments types
  
  // Data counts
  int<lower=0> N_pt;                      // total number of patients
  int<lower=0> N_res;                     // total number of responses
  
  // biomarker (interacting) covariant indicators
  matrix[N_pt, N_bm] bm_indicators;
  // treatment indicators
  int<lower=0, upper=1> tx_indicators[N_pt, N_tx];
  
  // response vectors
  vector<lower=0>[N_res] t_res;              // time (from diagnosis) of response
  int<lower=1, upper=N_pt> pt_index[N_res];  // index into patients
  vector[N_res] y_res;                       // response

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
  // total number of linear predictors
  int<lower=0> N_pred = N_int + N_bm + N_tx;
  // kludge to do math on indicators while enforcing binary data inputs
  matrix[N_pt, N_bm] bm_ind;
  matrix<lower=0, upper=1>[N_pt, N_tx] tx_ind;
  matrix[N_pt, N_int] int_ind;

  for (pt in 1:N_pt) {
    for (bm in 1:N_bm) {
      bm_ind[pt, bm] = bm_indicators[pt, bm];
    }
    for (tx in 1:N_tx) {
      tx_ind[pt, tx] = tx_indicators[pt, tx];
    }
    for (bm in 1:N_bm) {
      for (tx in 1:N_tx) {
	// [bm1-tx1, bm1-tx2, bm1-tx3, bm1-tx4, bm1-tx5, bm2-tx1, ..., bm5-tx5]
	int_ind[pt, (bm - 1) * N_tx + tx] = bm_indicators[pt, bm] * tx_indicators[pt, tx];
      }	
    }
  }
  
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

  vector[N_res] mu_res;
  
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
    
  //----------------------------------------------------------------------------
  // response sampling
  //----------------------------------------------------------------------------
  mu_res = mean_res(N_res, N_bm, N_tx, N_int,
		    t_res, pt_index, bm_ind, tx_ind, int_ind, 
		    beta0, beta_bm, beta_tx, beta_int);
  y_res ~ normal(mu_res, sigma_eps);
  
}
