/*
Multivariate Multi-level Generalized Linear Modeling of Cancer Patient Outcomes

This Stan program samples from the posterior probability distribution of the 
MM model, conditioned on time-series 
*/


functions {
  // project standard normal draws to a multivariate normal with zero mean
  matrix project_u(matrix z_u,      // standard normal draws, size (N_pt, N_re)
		   matrix L_u,      // cholesky decomposition, size (N_re, N_re)
		   vector sigma_u)  // diagonal scale terms, size N_re
  {
    int N_re = cols(z_u);
    int N_pt = rows(z_u);
    matrix[N_pt, N_re] u;
    u = (diag_pre_multiply(sigma_u, L_u) * z_u')';
    return u;
  }

  // inverse CDF sampling of HalfCauchy distribution
  real invcdf_hcauchy(real y,      // restricted to (0, 1), CDF value of HalfCauchy distribution
		      real scale)  // scale parameter of HalfCauchy distribution
  {
    return scale * tan(pi() * y / 2);
  }

  // broadcast a N_pt-sized vector up to size of the data vector
  vector rep_data_vector(vector x_pt,   // vector of size N_pt
			 int[] pt_idx)  // array of patient indices of size N_data
  {
    int N_pt = rows(x_pt);
    int N_data = num_elements(pt_idx);
    vector[N_data] x_data;
    for (i in 1:N_data) {
      x_data[i] = x_pt[pt_idx[i]];
    }
    return x_data;
  }

  // broadcast a (N_pt, N_re)-sized matrix up to size of the data matrix
  matrix rep_data_matrix(matrix x_pt,   // matrix of size (N_pt, N_re)
			 int[] pt_idx)  // array of patient indices of size N_data
  {
    int N_pt = rows(x_pt);
    int N_re = cols(x_pt);
    int N_data = num_elements(pt_idx);
    matrix[N_data, N_re] x_data;
    for (i in 1:N_data) {
      x_data[i] = x_pt[pt_idx[i]];
    }
    return x_data;
  }
  
  // population-level linear response
  vector dot_xbeta(matrix x,     // fixed-effect predictors, size (N_data, N_pred)
		   vector beta)	 // fixed-effect coefficients, size (N_pred,)     
  {
    int N_pred = rows(beta);
    int N_data = rows(x);
    matrix[N_data, N_pred] beta_data = rep_matrix(beta', N_data);
    return rows_dot_product(x, beta_data);
  }

  // patient-level linear response
  vector dot_zu(matrix z,      // random-effects predictors, size (N_data, N_re)
		matrix u,      // random-effects coefficients, size (N_pt, N_re)
		int[] pt_idx)  // index into patients, size (N_data,)
  {
    int N_data = rows(z);
    int N_re = cols(u);
    matrix[N_data, N_re] u_data = rep_data_matrix(u, pt_idx);
    return rows_dot_product(z, u_data);
  }
}


data {

  int N_pt;    // number of patients
  int N_bm;    // number of biomarkers
  int N_tx;    // number of treatments
  int N_tl;    // number of tumor load time series data
  
  matrix[N_pt, N_bm] bm_inputs;                     // biomarker inputs
  int<lower=0, upper=1> tx_indicators[N_pt, N_tx];  // treatment indicator array
  
  // TL data
  vector<lower=0>[N_tl] t_tl;                // time (from diagnosis) of TL measurement  
  int<lower=1, upper=N_pt> pt_idx_tl[N_tl];  // index into patients for TL measurements
  vector[N_tl] y_tl;                         // TL measurements

  // ###############################################################################################
  // Hyperparameters
  // ###############################################################################################

  // TL hyperparameters
  real loc_0_tl;
  real<lower=0> scale_0_tl;
  real loc_1_tl;
  real<lower=0> scale_1_tl;
  vector[N_bm] loc_bm_tl;
  vector<lower=0>[N_bm] scale_bm_tl;
  vector[N_tx] loc_tx_tl;
  vector<lower=0>[N_tx] scale_tx_tl;
  vector[N_bm * N_tx] loc_int_tl;
  vector<lower=0>[N_bm * N_tx] scale_int_tl;
  real<lower=0> scale_eps_tl;
  
  // Patient-level hyperparameters
  real<lower=0, upper=2> eta_u;
  vector<lower=0>[2] scale_u;

}


transformed data {

  int N_re = num_elements(scale_u);  // number of patient-level random effects
  int N_int = N_bm * N_tx;           // number of treatment-biomarker interactions

  // transformed data vectors
  vector[N_tl] log_y_tl = log(y_tl);
  
  // patient-level design matrices
  matrix[N_tl, 2] z_tl = append_col(rep_vector(1.0, N_tl), t_tl);

  // population-level design matrices
  // TL slope: x_tl = (1, x_bm, x_tx, x_int)
  matrix[N_pt, 1 + N_bm + N_tx + N_int] x_tl_pt;  // indexed by patient
  
  // kludge to do math on indicators while enforcing binary data inputs
  matrix<lower=0, upper=1>[N_pt, N_tx] tx_inputs;
  matrix[N_pt, N_int] int_inputs;
  for (pt in 1:N_pt) {
    for (tx in 1:N_tx) {
      tx_inputs[pt, tx] = tx_indicators[pt, tx];
      for (bm in 1:N_bm) {
	int_inputs[pt, (bm - 1) * N_tx + tx] = bm_inputs[pt, bm] * tx_inputs[pt, tx];
      }
    }
  }
  
  x_tl_pt = append_col(rep_vector(1.0, N_pt),
		       append_col(bm_inputs,
				  append_col(tx_inputs, int_inputs)));
}
  

parameters {
  
  // TL sub-model parameters
  real beta_0_tl;
  real beta_1_tl;  
  vector[N_bm] beta_bm_tl;
  vector[N_tx] beta_tx_tl;
  vector[N_bm * N_tx] beta_int_tl;
  real<lower=0, upper=1> s_eps_tl;

  // patient-level effects
  cholesky_factor_corr[N_re] L_u;
  vector<lower=0, upper=1>[N_re] s_u;
  matrix[N_pt, N_re] z_u;
  
}

transformed parameters {
  
  real sigma_eps_tl;
  vector[N_re] sigma_u;
  matrix[N_pt, N_re] u;

  vector[N_pt] delta_tl;
  vector[N_tl] eta_tl;
  
  sigma_eps_tl = invcdf_hcauchy(s_eps_tl, scale_eps_tl);

  for (i in 1:N_re) {
    sigma_u[i] = invcdf_hcauchy(s_u[i], scale_u[i]);
  }
  u = project_u(z_u, L_u, sigma_u);

  {
    vector[1 + N_bm + N_tx + N_int] beta_tl;
    vector[N_pt] u_tl_0 = u[1:N_pt, 1];
    vector[N_pt] u_tl_1 = u[1:N_pt, 2];
    
    // construct linear response for TL
    beta_tl = append_row(beta_1_tl,
			 append_row(beta_bm_tl,
				    append_row(beta_tx_tl, beta_int_tl)));
    delta_tl = dot_xbeta(x_tl_pt, beta_tl) + u_tl_1;
    eta_tl = rep_vector(beta_0_tl, N_tl)
      + rep_data_vector(u_tl_0, pt_idx_tl)
      + rep_data_vector(delta_tl, pt_idx_tl) .* t_tl;
  }
}


model {

  // ###############################################################################################
  // Priors
  // ###############################################################################################
  // 
  //   Noise parameters (sigma_eps_*, sigma_u) are implicitly drawn from a half-cauchy distribution
  //   by transforming a uniformly distributed parameter with the half-cauchy inverse CDF.
  // 
  //   Weibull shape parameters have a log-normal prior distribution by setting a normal prior over
  //   log(alpha).
  
  // TL population-level effects
  beta_0_tl ~ normal(loc_0_tl, scale_0_tl);
  beta_1_tl ~ normal(loc_1_tl, scale_1_tl);
  beta_bm_tl ~ normal(loc_bm_tl, scale_bm_tl);
  beta_tx_tl ~ normal(loc_tx_tl, scale_tx_tl);
  beta_int_tl ~ normal(loc_int_tl, scale_int_tl);

  // sample patient-level effects
  to_vector(z_u) ~ std_normal();

  // ###############################################################################################
  // Likelihoods
  // ###############################################################################################

  // likelihood for TL model
  log_y_tl ~ normal(eta_tl, rep_vector(sigma_eps_tl, N_tl));
  
}

generated quantities {
  
  corr_matrix[N_re] Omega_u;
  vector[N_tl] log_y_tl_ppc;
  
  Omega_u = L_u * L_u';
  log_y_tl_ppc = to_vector(normal_rng(eta_tl, rep_vector(sigma_eps_tl, N_tl)));

}
