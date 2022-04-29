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
  int N_ps;    // number of performance score time series data
  int N_tl;    // number of tumor load time series data
  
  matrix[N_pt, N_bm] bm_inputs;                     // biomarker inputs
  int<lower=0, upper=1> tx_indicators[N_pt, N_tx];  // treatment indicator array
  vector[N_pt] x_loc;                               // tumor location impact indicator
  
  // TL data
  vector<lower=0>[N_tl] t_tl;                // time (from diagnosis) of TL measurement  
  int<lower=1, upper=N_pt> pt_idx_tl[N_tl];  // index into patients for TL measurements
  vector[N_tl] y_tl;                         // TL measurements

  // PS data
  vector<lower=0>[N_ps] t_ps;                // time (from diagnosis) of PS measurement
  int<lower=1, upper=N_pt> pt_idx_ps[N_ps];  // index into patients for PS measurments
  vector[N_ps] y_ps;                         // PS measurements

  // SAE data
  vector[N_pt] T_sae;      // time-after-treatment of SAE
  int censored_sae[N_pt];  // 1 = censored, 0 = observed

  // DP data
  vector[N_pt] T_dp;       // time-after-treatment of DP
  int censored_dp[N_pt];   // 1 = censored, 0 = observed

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
  
  // DP hyperparameters
  real loc_gamma_tl_dp;
  real<lower=0> scale_gamma_tl_dp;
  real loc_log_alpha_dp;
  real<lower=0> scale_log_alpha_dp;

  // SAE hyperparameters
  vector[N_tx] loc_tx_sae;
  vector<lower=0>[N_tx] scale_tx_sae;
  vector[N_bm * N_tx] loc_int_sae;
  vector<lower=0>[N_bm * N_tx] scale_int_sae;
  real loc_log_alpha_sae;
  real<lower=0> scale_log_alpha_sae;

  // PS hyperparameters
  real loc_0_ps;
  real<lower=0> scale_0_ps;
  real loc_1_ps;
  real<lower=0> scale_1_ps;
  vector[N_tx] loc_tx_ps;
  vector<lower=0>[N_tx] scale_tx_ps;
  vector[N_tx] loc_sae_ps;
  vector<lower=0>[N_tx] scale_sae_ps;  
  real<lower=0> scale_eps_ps;

  // Patient-level hyperparameters
  real<lower=0, upper=2> eta_u;
  vector<lower=0>[5] scale_u;

}


transformed data {

  int N_re = num_elements(scale_u);  // number of patient-level random effects
  int N_int = N_bm * N_tx;           // number of treatment-biomarker interactions

  // transformed data vectors
  vector[N_tl] log_y_tl = log(y_tl);
  vector[N_ps] logit_y_ps = logit(y_ps);
  
  // patient-level design matrices
  matrix[N_tl, 2] z_tl = append_col(rep_vector(1.0, N_tl), t_tl);
  matrix[N_ps, 2] z_ps = append_col(rep_vector(1.0, N_ps), t_ps);

  // population-level design matrices
  // TL slope: x_tl = (1, x_bm, x_tx, x_int)
  matrix[N_pt, 1 + N_bm + N_tx + N_int] x_tl_pt;  // indexed by patient
  matrix[N_tl, 1 + N_bm + N_tx + N_int] x_tl_tl;  // indexed by t_tl
  matrix[N_ps, 1 + N_bm + N_tx + N_int] x_tl_ps;  // indexed by t_ps
  // SAE "slope": x_sae = (x_tx, x_int)
  matrix[N_pt, N_tx + N_int] x_sae_pt;            // indexed by patient
  // treatment indicators
  matrix[N_ps, N_tx] x_tx_ps;                     // indexed by t_ps
  // impactful location indicators
  vector[N_ps] x_loc_ps;                          // indexed by t_ps
  // SAE impact indicator
  matrix[N_ps, N_tx] I_sae;
  
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

  // SAE impact indicators
  for (j in 1:N_ps) {
    int pt = pt_idx_ps[j];
    for (tx in 1:N_tx) {
      if ((t_ps[j] > T_sae[pt]) && (tx_inputs[pt, tx]) && !(censored_sae[pt])) {
	I_sae[j, tx] = 1.0;
      }
      else {
	I_sae[j, tx] = 0.0;
      }
    }
  }
  
  x_tl_pt = append_col(rep_vector(1.0, N_pt),
		       append_col(bm_inputs,
				  append_col(tx_inputs, int_inputs)));
  x_tl_tl = rep_data_matrix(x_tl_pt, pt_idx_tl);
  x_tl_ps = rep_data_matrix(x_tl_pt, pt_idx_ps);
  x_sae_pt = append_col(tx_inputs, int_inputs);
  x_tx_ps = rep_data_matrix(tx_inputs, pt_idx_ps);
  x_loc_ps = rep_data_vector(x_loc, pt_idx_ps);
}
  

parameters {
  
  // TL sub-model parameters
  real beta_0_tl;
  real beta_1_tl;  
  vector[N_bm] beta_bm_tl;
  vector[N_tx] beta_tx_tl;
  vector[N_bm * N_tx] beta_int_tl;
  real<lower=0, upper=1> s_eps_tl;

  // DP sub-model parameters
  real gamma_tl_dp;
  real log_alpha_dp;
  
  // SAE sub-model parameters
  vector[N_tx] beta_tx_sae;
  vector[N_bm * N_tx] beta_int_sae;
  real log_alpha_sae;
  
  // PS sub-model population-effects
  real beta_0_ps;
  real beta_1_ps;
  vector[N_tx] beta_tx_ps;
  real beta_tl_ps;
  vector[N_tx] beta_sae_ps;
  real<lower=0, upper=1> s_eps_ps;

  // patient-level effects
  cholesky_factor_corr[N_re] L_u;
  vector<lower=0, upper=1>[N_re] s_u;
  matrix[N_pt, N_re] z_u;
  
}

transformed parameters {
  
  real sigma_eps_tl;
  real sigma_eps_ps;
  real alpha_dp;
  real alpha_sae;
  vector[N_re] sigma_u;
  matrix[N_pt, N_re] u;
  
  sigma_eps_tl = invcdf_hcauchy(s_eps_tl, scale_eps_tl);
  sigma_eps_ps = invcdf_hcauchy(s_eps_ps, scale_eps_ps);
  alpha_dp = exp(log_alpha_dp);
  alpha_sae = exp(log_alpha_sae);

  for (i in 1:N_re) {
    sigma_u[i] = invcdf_hcauchy(s_u[i], scale_u[i]);
  }
  u = project_u(z_u, L_u, sigma_u);
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

  // DP parameters
  gamma_tl_dp ~ normal(loc_gamma_tl_dp, scale_gamma_tl_dp);
  log_alpha_dp ~ normal(loc_log_alpha_dp, scale_log_alpha_dp);
  
  // SAE population-level effects
  beta_tx_sae ~ normal(loc_tx_sae, scale_tx_sae);
  beta_int_sae ~ normal(loc_int_sae, scale_int_sae);
  // SAE weibull shape parameter
  log_alpha_sae ~ normal(loc_log_alpha_sae, scale_log_alpha_dp);
  
  // PS sub-model population-effects
  beta_0_ps ~ normal(loc_0_ps, scale_0_ps);
  beta_1_ps ~ normal(loc_1_ps, scale_1_ps);
  beta_tx_ps ~ normal(loc_tx_ps, scale_tx_ps);
  // beta_int_ps ~ normal(loc_int_ps, scale_int_ps);
  beta_sae_ps ~ normal(loc_sae_ps, scale_sae_ps);

  // sample patient-level effects
  to_vector(z_u) ~ std_normal();

  // ###############################################################################################
  // Likelihoods
  // ###############################################################################################

  {
    vector[1 + N_bm + N_tx + N_int] beta_tl;
    vector[N_tx + N_int] beta_sae;
    vector[3 + 2 * N_tx] beta_ps;

    real sigma_dp;
    real sigma_sae;

    matrix[N_ps, 3 + 2 * N_tx] x_ps;

    vector[N_pt] u_tl_0;
    vector[N_pt] u_tl_1;
    matrix[N_pt, 2] u_ps;
    vector[N_pt] u_sae;

    vector[N_pt] delta_tl;
    vector[N_pt] delta_sae;
    vector[N_tl] eta_tl;
    vector[N_ps] eta_tl_ps;
    vector[N_ps] eta_ps;

    u_tl_0 = u[1:N_pt, 1];
    u_tl_1 = u[1:N_pt, 2];
    u_ps = u[1:N_pt, 3:4];
    u_sae = u[1:N_pt, 5];

    // construct linear responses for TL, DP and SAE
    beta_tl = append_row(beta_1_tl,
			 append_row(beta_bm_tl,
				    append_row(beta_tx_tl, beta_int_tl)));
    beta_sae = append_row(beta_tx_sae, beta_int_sae);
    delta_tl = dot_xbeta(x_tl_pt, beta_tl) + u_tl_1;
    eta_tl = rep_vector(beta_0_tl, N_tl)
             + rep_data_vector(u_tl_0, pt_idx_tl)
             + rep_data_vector(delta_tl, pt_idx_tl) .* t_tl;
    delta_sae = dot_xbeta(x_sae_pt, beta_sae) + u_sae;

    // construct linear responses for PS
    beta_ps = append_row(beta_0_ps,
			 append_row(beta_tx_ps,
				    append_row(beta_tl_ps,
					       append_row(beta_sae_ps, beta_1_ps))));
    eta_tl_ps = rep_vector(beta_0_tl, N_ps)
                + rep_data_vector(u_tl_0, pt_idx_ps)
                + rep_data_vector(delta_tl, pt_idx_ps) .* t_ps;
    x_ps = append_col(rep_vector(1.0, N_ps),
		      append_col(x_tx_ps,
				 append_col(x_loc_ps .* eta_tl_ps,
					    append_col(I_sae, t_ps))));
    eta_ps = dot_xbeta(x_ps, beta_ps) + dot_zu(z_ps, u_ps, pt_idx_ps);

    // likelihood for TL model
    log_y_tl ~ normal(eta_tl, rep_vector(sigma_eps_tl, N_tl));

    // likelihood for PS model
    logit_y_ps ~ normal(eta_ps, rep_vector(sigma_eps_ps, N_ps));

    // likelihoods for DP and SAE models
    for (pt in 1:N_pt) {
      // DP model
      sigma_dp = exp(gamma_tl_dp * delta_tl[pt]) ^ (-1 / alpha_dp);
      if (censored_dp[pt]) {
	target += weibull_lccdf(T_dp[pt] | alpha_dp, sigma_dp);
      }
      else {
	T_dp[pt] ~ weibull(alpha_dp, sigma_dp);
      }
      // SAE model
      sigma_sae = exp(delta_sae[pt]) ^ (-1 / alpha_sae);
      if (censored_sae[pt]) {
	target += weibull_lccdf(T_sae[pt] | alpha_sae, sigma_sae);
      }
      else {
	T_sae[pt] ~ weibull(alpha_sae, sigma_sae);
      }
    }
  }
}

generated quantities {
  
  corr_matrix[N_re] Omega_u;
  Omega_u = L_u * L_u';

}
