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
		  vector beta_int, // interaction effects term
		  matrix u)        // N_pt by 2 matrix of random effects intercepts and slopes
  {
    vector[N_res] intercept;
    vector[N_res] random_effects_slope;
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
    
      intercept[i] = beta0 + u[pt, 1];
      random_effects_slope[i] = u[pt, 2] * t;

      // build predictors
      x_bm = (bm_ind[pt, :] * t)';
      x_tx = (tx_ind[pt, :] * t)';
      x_int = (int_ind[pt, :] * t)';
    
      fixed_effects_slope[i] =
	+ dot_product(beta_bm, x_bm)
	+ dot_product(beta_tx, x_tx)
	+ dot_product(beta_int, x_int);
    }
    return intercept + random_effects_slope + fixed_effects_slope;
  }


  matrix project_u(matrix z_u,    // standard normal draws
		   real mu_u1,    // mean slope
		   vector tau_u,  // scale vector
		   real rho_u)    // correlation
  {
    vector[2] mu_u = [0.0, mu_u1]';
    int N_pt = rows(z_u);
    matrix[2, 2] L_u = [[1.0, 0.0], [rho_u, sqrt(1.0 - rho_u ^ 2)]];
    matrix[N_pt, 2] u = rep_matrix(mu_u', N_pt) + (diag_pre_multiply(tau_u, L_u) * z_u')';
    return u;
  }

}


data {

  int<lower=0> N_samples; // number of parameter space samples
  
  // Parameter counts
  int<lower=0> N_bm;  // number of biomarkers (with treatment interaction effects)
  int<lower=0> N_tx;  // number of treatments types
  
  // Data counts
  int<lower=0> N_pt;                      // total number of patients
  int<lower=0> N_res;                     // total number of responses
  
  // biomarker (interacting) covariant indicators
  int<lower=0, upper=1> bm_indicators[N_pt, N_bm];
  // treatment indicators
  int<lower=0, upper=1> tx_indicators[N_pt, N_tx];
  
  // response vectors
  vector<lower=0>[N_res] t_res;              // time (from diagnosis) of response
  int<lower=1, upper=N_pt> pt_index[N_res];  // index into patients

  // hyperparameters
  
  // parameters
  // fixed (biomarker/tx-level) effects
  vector<lower=0, upper=100>[N_samples] beta0;  // intercept
  matrix[N_samples, N_bm] beta_bm;              // biomarker effects
  matrix[N_samples, N_tx] beta_tx;              // treatment effects
  matrix[N_samples, N_bm * N_tx] beta_int;      // interaction effects

  // random (patient-level) effects hyperparameters
  vector[N_samples] mu_u1;                    // mean of intercept and slope
  matrix<lower=0>[N_samples, 2] tau_u;        // scale parameters for intercept and slope
  vector<lower=-1, upper=1>[N_samples] rho_u; // correlation of intercept and slope
    
  // measurement standard deviation  
  vector<lower=0>[N_samples] sigma_eps;
  
}


transformed data {
  // number of treatment-biomarker interactions
  int<lower=0> N_int = N_bm * N_tx;  
  // total number of linear predictors
  int<lower=0> N_pred = N_int + N_bm + N_tx;
  // kludge to do math on indicators while enforcing binary data inputs
  matrix<lower=0, upper=1>[N_pt, N_bm] bm_ind;
  matrix<lower=0, upper=1>[N_pt, N_tx] tx_ind;
  matrix<lower=0, upper=1>[N_pt, N_int] int_ind;  
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
  matrix[N_pt, 2] z_u; // uncorrelated random effects
}


model {
  to_vector(z_u) ~ std_normal();
}


generated quantities {
  
  vector[N_res] y_res;  // predicted response

  {
    vector[N_res] mu_res; // mean response
    int j;                // index into samples
    matrix[N_pt, 2] u;    // patient level effects
    
    j = categorical_rng(rep_vector(inv(N_samples), N_samples));
    u = project_u(z_u, mu_u1[j]', tau_u[j]', rho_u[j]);
    
    mu_res = mean_res(N_res, N_bm, N_tx, N_int, t_res, pt_index,
		      bm_ind, tx_ind, int_ind,
		      beta0[j]', beta_bm[j]', beta_tx[j]', beta_int[j]', u);
    y_res = to_vector(normal_rng(mu_res, sigma_eps[j]));
  }
}
