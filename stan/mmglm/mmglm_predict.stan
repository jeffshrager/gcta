/*
Multivariate Multi-level Generalized Linear Modeling of Cancer Patient Outcomes

This Stan program generates mock data by sampling from the predictive
distribution of the MM model, marginalized over the posterior probability
distribution.
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

  // inverse CDF sampling of Cauchy distribution
  real invcdf_cauchy(real y,      // restricted to (0, 1), CDF value of Cauchy distribution
		     real loc,    // location parameter of Cauchy distribution
		     real scale)  // scale parameter of Cauchy distribution
  {
    return loc + scale * tan(pi() * (y - 0.5));
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
  vector<lower=0>[N_ps] t_tl;                // time (from diagnosis) of TL measurement  
  int<lower=1, upper=N_pt> pt_idx_tl[N_tl];   // index into patients for TL measurements

  // PS data
  vector<lower=0>[N_ps] t_ps;               // time (from diagnosis) of PS measurement
  int<lower=1, upper=N_pt> pt_idx_ps[N_ps]; // index into patients for PS measurments
  
  // parameters
  int N_samples;  // number of parameter samples

  // TL sub-model parameters
  vector[N_samples] beta_0_tl;
  vector[N_samples] beta_1_tl;  
  matrix[N_samples, N_bm] beta_bm_tl;
  matrix[N_samples, N_tx] beta_tx_tl;
  matrix[N_samples, N_bm * N_tx] beta_int_tl;
  vector[N_samples] sigma_eps_tl;

  // DP sub-model parameters
  vector[N_samples] gamma_tl_dp;
  vector[N_samples] log_alpha_dp;
  
  // SAE sub-model parameters
  matrix[N_samples, N_tx] beta_tx_sae;
  matrix[N_samples, N_bm * N_tx] beta_int_sae;
  vector[N_samples] log_alpha_sae;
  
  // PS sub-model population-effects
  vector[N_samples] beta_0_ps;
  vector[N_samples] beta_1_ps;
  matrix[N_samples, N_tx] beta_tx_ps;
  vector[N_samples] beta_tl_ps;
  matrix[N_samples, N_tx] beta_sae_ps;
  vector[N_samples] sigma_eps_ps;

  // patient-level effects
  corr_matrix[5] Omega_u[N_samples];
  matrix<lower=0>[N_samples, 5] sigma_u;
}


transformed data {

  int N_re = cols(sigma_u);                    // number of patient-level random effects
  int N_int = N_bm * N_tx;                     // number of treatment-biomarker interactions
  
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
  x_tl_tl = rep_data_matrix(x_tl_pt, pt_idx_tl);
  x_tl_ps = rep_data_matrix(x_tl_pt, pt_idx_ps);
  x_sae_pt = append_col(tx_inputs, int_inputs);
  x_tx_ps = rep_data_matrix(tx_inputs, pt_idx_ps);
  x_loc_ps = rep_data_vector(x_loc, pt_idx_ps);
}
  

parameters {
  matrix[N_pt, N_re] z_u;
}


model {
  to_vector(z_u) ~ std_normal();
}

generated quantities {

  vector[N_tl] y_tl;
  vector[N_ps] y_ps;
  vector[N_pt] T_sae;
  vector[N_pt] T_dp;
  {
    int i;                  // index into parameter samples
    int pt;                 // index into patients
    matrix[N_pt, N_re] u;   // patient-level effects
    matrix[N_re, N_re] L_u; // cholesky-decomposition of correlation matrix
    
    vector[1 + N_bm + N_tx + N_int] beta_tl;
    vector[N_tx + N_int] beta_sae;
    vector[3 + 2 * N_tx] beta_ps;

    real alpha;
    real sigma;

    matrix[N_ps, N_tx] I_sae;
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
    vector[N_ps] logit_y_ps;
    
    // sample from parameter distribution
    i = categorical_rng(rep_vector(inv(N_samples), N_samples));
    
    // patient-level effects
    L_u = cholesky_decompose(Omega_u[i]);
    u = project_u(z_u, L_u, sigma_u[i]');
    u_tl_0 = u[1:N_pt, 1];
    u_tl_1 = u[1:N_pt, 2];    
    u_ps = u[1:N_pt, 3:4];
    u_sae = u[1:N_pt, 5];
	   
    // construct linear responses for TL and SAE
    beta_tl = append_row(beta_1_tl[i]',
			 append_row(beta_bm_tl[i]',
				    append_row(beta_tx_tl[i]', beta_int_tl[i]')));
    beta_sae = append_row(beta_tx_sae[i]', beta_int_sae[i]');
    delta_tl = dot_xbeta(x_tl_pt, beta_tl) + u_tl_1;
    eta_tl = rep_vector(beta_0_tl[i], N_tl)
             + rep_data_vector(u_tl_0, pt_idx_tl)
             + rep_data_vector(delta_tl, pt_idx_tl) .* t_tl;
    delta_sae = dot_xbeta(x_sae_pt, beta_sae) + u_sae;

    // sample from TL model
    y_tl = to_vector(lognormal_rng(eta_tl, rep_vector(sigma_eps_tl[i], N_tl)));
    
    // sample from time-to-disease progression model
    alpha = exp(log_alpha_dp[i]);
    for (j in 1:N_pt) {
      sigma = exp(gamma_tl_dp[i] * delta_tl[j]) ^ (-1 / alpha);
      T_dp[j] = weibull_rng(alpha, sigma);
    }
    
    // sample from time-to-SAE model
    alpha = exp(log_alpha_sae[i]);
    for (j in 1:N_pt) {
      sigma = exp(delta_sae[j]) ^ (-1 / alpha);
      T_sae[j] = weibull_rng(alpha, sigma);
    }

    // construct linear responses for PS
    beta_ps = append_row(beta_0_ps[i],
			 append_row(beta_tx_ps[i, 1:N_tx]',
				    append_row(beta_tl_ps[i],
					       append_row(beta_sae_ps[i, 1:N_tx]', beta_1_ps[i]))));
    // SAE impact indicators
    for (j in 1:N_ps) {
      pt = pt_idx_ps[j];
      for (tx in 1:N_tx)
	{
	  if ((t_ps[j] > T_sae[pt]) && (tx_inputs[pt, tx])) {
	    I_sae[j, tx] = 1.0;
	  }
	  else {
	    I_sae[j, tx] = 0.0;
	  }
	}
    }

    eta_tl_ps = rep_vector(beta_0_tl[i], N_ps)
                + rep_data_vector(u_tl_0, pt_idx_ps)
                + rep_data_vector(delta_tl, pt_idx_ps) .* t_ps;
    x_ps = append_col(rep_vector(1.0, N_ps),
		      append_col(x_tx_ps,
				 append_col(x_loc_ps .* eta_tl_ps,
					    append_col(I_sae, t_ps))));
    eta_ps = dot_xbeta(x_ps, beta_ps) + dot_zu(z_ps, u_ps, pt_idx_ps);

    // sample from PS model
    logit_y_ps = to_vector(normal_rng(eta_ps, rep_vector(sigma_eps_ps[i], N_ps)));
    y_ps = inv_logit(logit_y_ps);
  }
}


