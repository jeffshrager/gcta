functions {
  // project standard normal draws to a multivariate normal with zero mean
  matrix project_u(matrix z_u,       // standard normal draws, size (N_pt, N_re)
		   matrix L_u,       // cholesky decomposition, size (N_re, N_re)
		   vector sigma_u)   // diagonal scale terms, size N_re
  {
    int N_re = rows(L_u);
    int N_pt = rows(z_u);
    matrix[N_pt, N_re] u;
    u = (diag_pre_multiply(sigma_u, L_u) * z_u')';
    return u;
  }

  // inverse CDF sampling of Cauchy distribution
  real invcdf_cauchy(real y,      // restricted to (0, 1), CDF value of Cauchy distribution
		     real loc,    // location parameter of Cauchy distribution
		     real scale)  // scale parameter of Cauchy distribution
  {
    return loc + scale * tan(pi() * (y - 0.5));
  }

  // mean linear response for the population- and patient-level predictors/coefficients
  vector compute_eta(matrix x,     // fixed-effect predictors, size (N_data, N_pred)
		     vector beta,  // fixed-effect coefficients, size (N_pred,)
		     matrix z,     // random-effects predictors, size (N_data, 2)
		     matrix u,     // random-effects coefficients, size (N_pt, 2)
		     int[] pt_idx) // index into patients, size (N_data,)
  {
    int N_pred = num_elements(beta);
    int N_data = rows(x);
    int N_pt = cols(u);
    
    int pt;
    matrix[N_data, N_pred] beta_data;
    matrix[N_data, 2] u_data;
    vector[N_data] eta;

    beta_data = rep_matrix(beta', N_data);
    for (i in 1:N_data) {
      pt = pt_idx[i];
      u_data[i, 1:2] = u[pt, 1:2];
    }
    eta = rows_dot_product(x, beta_data) + rows_dot_product(z, u_data);
    return eta;
  }

  // hazard function for Weibull surivival model
  vector lambda_weibull(vector T,    // time
			vector eta,  // log hazard ratio
			real k)      // shape parameter
  {
    int N = num_elements(T);
    vector[N] lambda;
    for (i in 1:N) {
      lambda[i] = k * (T[i] ^ (k - 1)) * exp(-eta[i]);
    }
    return lambda;
  }
  // survival function for Weibull survival model
  vector S_weibull(vector T,    // time
		   vector eta,  // log hazard ratio
		   real k)      // shape parameter
  {
    int N = num_elements(T);
    vector[N] S;
    for (i in 1:N) {
      S[i] = exp(-eta[i]) * (1 - exp(-(T[i] ^ k)));
    }
    return S;
  }
  // log probability of right-censored Weibull survival times
  vector censored_weibull_lp(vector T_obs,    // observed survival times, size (N_pt,)
			     int[] censored,  // censored (1) or uncensored (0), size (N_pt,)
			     vector eta,      // log hazard ratio, size (N_pt,)
			     real k)          // shape parameter

  {
    int N_pt = num_elements(T_obs);
    vector[N_pt] lp;

    vector[N_pt] S = S_weibull(T_obs, eta, k);
    vector[N_pt] lambda = lambda_weibull(T_obs, eta, k);

    for (pt in 1:N_pt) {
      if(censored[pt]) {
	lp[pt] = log(1 - S[pt]);
      }
      else {
	lp[pt] = log(S[pt]) + log(lambda[pt]);
      }
    }
    return lp;
  }
  
}

data {
  int N_pt;    // number of patients
  int N_bm;    // number of biomarkers
  int N_tx;    // number of treatments
  int N_kps;   // number of Karnofsky Performance score time-series data
  int N_tl;    // number of Tumor Load time-series data
  
  matrix[N_pt, N_bm] bm_inputs;                     // biomarker inputs
  int<lower=0, upper=1> tx_indicators[N_pt, N_tx];  // treatment indicator array
  
  // TL data
  vector<lower=0>[N_kps] t_tl;                // time (from diagnosis) of TL measurement  
  int<lower=1, upper=N_pt> pt_idx_tl[N_tl];   // index into patients for TL measurements
  vector[N_tl] y_tl;                          // TL measurements

  // KPS data
  vector<lower=0>[N_kps] t_kps;               // time (from diagnosis) of KPS measurement
  int<lower=1, upper=N_pt> pt_idx_kps[N_kps]; // index into patients for KPS measurments
  vector[N_kps] y_kps;                        // KPS measurements

  // Survival data
  vector[N_pt] T_obs;                         // Observed survival time
  int<lower=0, upper=1> censored[N_pt];       // censored (1) or non-censored (0)

  // patient-level effects hyperparameters
  real eta_u;
  vector<lower=0>[6] scale_sigma_u;

  // population-level effects hyperparameters
  // for TL
  vector[2] loc_const_tl;
  cov_matrix[2] scale_const_tl;
  vector[N_bm] loc_bm_tl;
  cov_matrix[N_bm] scale_bm_tl;
  vector[N_tx] loc_tx_tl;
  cov_matrix[N_tx] scale_tx_tl;
  vector[N_bm * N_tx] loc_int_tl;
  cov_matrix[N_bm * N_tx] scale_int_tl;
  real<lower=0> scale_sigma_eps_tl;
  
  // for KPS
  vector[2] loc_const_kps;
  cov_matrix[2] scale_const_kps;
  real loc_tl_kps;
  real scale_tl_kps;
  vector[N_bm] loc_bm_kps;
  cov_matrix[N_bm] scale_bm_kps;
  vector[N_tx] loc_tx_kps;
  cov_matrix[N_tx] scale_tx_kps;
  vector[N_bm * N_tx] loc_int_kps;
  cov_matrix[N_bm * N_tx] scale_int_kps;
  real<lower=0> scale_sigma_eps_kps;
  
  // for survival model
  vector[2] loc_const_surv;
  cov_matrix[2] scale_const_surv;
  real loc_tl_surv;
  real scale_tl_surv;
  real loc_kps_surv;
  real scale_kps_surv;
  real<lower=0> scale_sigma_eps_surv;
}


transformed data {

  int N_re = num_elements(scale_sigma_u);  // number of patient-level random effects
  int N_int = N_bm * N_tx;                 // number of treatment-biomarker interactions
  int N_pred = 2 + N_bm + N_tx + N_int;    // number of base predictors
  real kps_max = 100.0;                    // maximum of KPS
  // log transform of TL data
  vector[N_tl] log_y_tl = log(y_tl);
  
  // patient-level design matrices
  matrix[N_tl, 2] z_tl = append_col(rep_vector(1.0, N_tl), t_tl);
  matrix[N_kps, 2] z_kps = append_col(rep_vector(1.0, N_kps), t_kps);
  matrix[N_pt, 2] z_surv = append_col(rep_vector(1.0, N_pt), T_obs);

  // population-level design matrices
  matrix[N_tl, N_pred] x_tl;
  matrix[N_kps, N_pred] x_kps;
  
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

  // build population-level design matrices
  {
    int start_col;
    int end_col;
    real t;
    int pt;
    // predictors for TL
    for (i in 1:N_tl) {
      pt = pt_idx_tl[i];
      t = t_tl[i];
      start_col = 1;
      x_tl[i, 1] = 1.0;
      x_tl[i, 2] = t;
      start_col = 3;
      end_col = start_col + N_bm - 1;
      x_tl[i, start_col:end_col] = bm_inputs[pt, 1:N_bm] * t;
      start_col = end_col;
      end_col = start_col + N_tx - 1;
      x_tl[i, start_col:end_col] = tx_inputs[pt, 1:N_tx] * t;
      start_col = end_col;
      end_col = start_col + N_int - 1;
      x_tl[i, start_col:end_col] = int_inputs[pt, 1:N_int] * t;
    }
    // predictors for KPS
    for (i in 1:N_kps) {
      pt = pt_idx_kps[i];
      t = t_kps[i];
      start_col = 1;
      x_kps[i, 1] = 1.0;
      x_kps[i, 2] = t;
      start_col = 3;
      end_col = start_col + N_bm - 1;
      x_kps[i, start_col:end_col] = bm_inputs[pt, 1:N_bm] * t;
      start_col = end_col;
      end_col = start_col + N_tx - 1;
      x_kps[i, start_col:end_col] = tx_inputs[pt, 1:N_tx] * t;
      start_col = end_col;
      end_col = start_col + N_int - 1;
      x_kps[i, start_col:end_col] = int_inputs[pt, 1:N_int] * t;
    }
  }  
}


parameters {
  // TL sub-model population-effects
  vector[2] beta_const_tl;
  vector[N_bm] beta_bm_tl;
  vector[N_tx] beta_tx_tl;
  vector[N_int] beta_int_tl;
  real<lower=0, upper=1> s_tl;

  // KPS sub-model population-effects
  vector[2] beta_const_kps;
  vector[N_bm] beta_bm_kps;
  vector[N_tx] beta_tx_kps;
  vector[N_int] beta_int_kps;
  real beta_tl_kps;
  real<lower=0, upper=1> s_kps;

  // Survival sub-model
  real<lower=-10, upper=10> log_k_weibull;
  vector[2] beta_const_surv;
  real beta_tl_surv;
  real beta_kps_surv;
  real<lower=0, upper=1> s_surv;

  // patient-level effects
  matrix[N_pt, N_re] z_u;
  cholesky_factor_corr[N_re] L_u;
  vector<lower=0, upper=1>[N_re] s_u;
}


transformed parameters {
  real sigma_eps_tl;
  real sigma_eps_kps;
  real sigma_eps_surv;
  real k_weibull;
  vector[N_re] sigma_u;
  matrix[N_pt, N_re] u;

  k_weibull = exp(log_k_weibull);
  for (i in 1:N_re) {
    sigma_u[i] = invcdf_cauchy(s_u[i], 0.0, scale_sigma_u[i]);
  }
  sigma_eps_tl = invcdf_cauchy(s_tl, 0.0, scale_sigma_eps_tl);
  sigma_eps_kps = invcdf_cauchy(s_kps, 0.0, scale_sigma_eps_kps);
  sigma_eps_surv = invcdf_cauchy(s_surv, 0.0, scale_sigma_eps_surv);
  u = project_u(z_u, L_u, sigma_u);
}


model {

  vector[N_pred] beta_tl;
  vector[N_pred] beta_kps;
  
  vector[N_tl] eta_tl;
  vector[N_kps] eta_kps;
  vector[N_pt] eta_surv;
  
  matrix[N_pt, 2] u_tl = u[1:N_pt, 1:2];
  matrix[N_pt, 2] u_kps = u[1:N_pt, 3:4];
  matrix[N_pt, 2] u_surv = u[1:N_pt, 5:6];
  
  // priors
  // population-level TL sub-model
  beta_const_tl ~ multi_normal(loc_const_tl, scale_const_tl);
  beta_bm_tl ~ multi_normal(loc_bm_tl, scale_bm_tl);
  beta_tx_tl ~ multi_normal(loc_tx_tl, scale_tx_tl);
  beta_int_tl ~ multi_normal(loc_int_tl, scale_int_tl);
  // population-level KPS sub-model
  beta_const_kps ~ multi_normal(loc_const_kps, scale_const_kps);
  beta_bm_kps ~ multi_normal(loc_bm_kps, scale_bm_kps);
  beta_tx_kps ~ multi_normal(loc_tx_kps, scale_tx_kps);
  beta_int_kps ~ multi_normal(loc_int_kps, scale_int_kps);
  beta_tl_kps ~ normal(loc_tl_kps, scale_tl_kps);
  // population-level survival sub-model
  beta_const_kps ~ multi_normal(loc_const_surv, scale_const_surv);
  beta_tl_surv ~ normal(loc_tl_surv, scale_tl_surv);
  beta_kps_surv ~ normal(loc_kps_surv, scale_kps_surv);
  // patient-level
  // sigma_u is implicitly drawn from a Cauchy distribution by inverse CDF sampling
  to_vector(z_u) ~ std_normal();
  L_u ~ lkj_corr_cholesky(eta_u);

  // construct linear responses
  beta_tl = append_row(beta_const_tl,
		       append_row(beta_bm_tl,
				  append_row(beta_tx_tl, beta_int_tl)));
  beta_kps = append_row(beta_const_kps,
		       append_row(beta_bm_kps,
				  append_row(beta_tx_kps, beta_int_kps)));
  eta_tl = compute_eta(x_tl, beta_tl, z_tl, u_tl, pt_idx_tl);
  eta_kps = compute_eta(x_kps, beta_kps, z_kps, u_kps, pt_idx_kps) + beta_tl_kps * eta_tl;
  eta_surv = beta_tl_surv * eta_tl + beta_kps_surv * eta_kps + rows_dot_product(z_surv, u_surv);

  // longitudinal likelihoods
  log_y_tl ~ normal(eta_tl, sigma_eps_tl);
  y_kps ~ logistic(kps_max * logistic_cdf(eta_kps, 0.0, 1.0), sigma_eps_kps);
  // random-effects likelihood is already implicit as a prior over z_u, projected to u
  // survival likelihood
  target += sum(censored_weibull_lp(T_obs, censored, eta_surv, k_weibull));
}


/* generated quantities { */

/* } */
