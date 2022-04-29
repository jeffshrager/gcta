functions {
  
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
  // input data
  int N_pt;                              // number of patients		 
  int N_data;                            // number of data points in total 
  int N_tx;                              // number of treatments		 
  int<lower=1, upper=N_pt> idx[N_data];  // patient index into data vectors
  matrix[N_data, 1 + N_tx] x;            // predictors for fixed-effects (1, t, tx_A, tx_B)
  matrix[N_data, 2] z;                   // predictors for random-effects (1, t)
  vector[N_data] y;                      // data vector

  // hyperparameters for prior
  vector[1 + N_tx] loc_beta;
  vector[1 + N_tx] scale_beta;
  real scale_eps;
  // real eta_u;
  vector[2] scale_u;
}


parameters {
  vector[1 + N_tx] beta;             // fixed-effect slopes
  real<lower=0.5, upper=1> s_eps;    // standardized noise scale
  vector<lower=0.5, upper=1>[2] s_u; // standardized random-effects scale
  matrix[N_pt, 2] z_u;               // standardized random-effects terms
}


transformed parameters {
  real sigma_eps;
  vector[2] sigma_u;
  matrix[N_pt, 2] u;
  {
    sigma_eps = invcdf_cauchy(s_eps, 0.0, scale_eps);
    for (i in 1:2) {
      sigma_u[i] = invcdf_cauchy(s_u[i], 0.0, scale_u[i]);
    }
    u = project_u(z_u, append_row([1, 0], [0, 1]), scale_u);
  }
}


model {
  vector[N_data] eta;
  // priors
  //   sigma_eps, sigma_u are implicitly drawn from a Cauchy prior
  //   via a parameter tranformation over a uniform distribution
  beta ~ normal(loc_beta, scale_beta);
  to_vector(z_u) ~ std_normal();

  // likelihood
  eta = dot_xbeta(x, beta) + dot_zu(z, u, idx);
  y ~ normal(eta, sigma_eps);
}

generated quantities {
  vector[N_data] y_ppc;
  {
    vector[N_data] eta = dot_xbeta(x, beta) + dot_zu(z, u, idx);
    for (i in 1:N_data) {
      y_ppc[i] = normal_rng(eta[i], sigma_eps);
    }
  }
}
