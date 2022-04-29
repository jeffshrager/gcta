import json

import numpy as np

randn = np.random.randn

with open("patients.py") as f:
    exec(f.read())
    
N_samples = 1

time_scale = 6    # months
tl_scale = 1e-2   # delta TL per month
dp_scale = 1e2   # delta DP hazard per TL per month
sae_scale = 1e1  # delta SAE hazard per month
ps_scale = 1e-2   # delta PS per month
eps_scale = 10
u_scale = 2 

################################################################################
# TL sub-model
################################################################################

beta_0_tl = [0]
beta_1_tl = [0]
beta_bm_tl = [(10 * tl_scale * np.ones(N_bm)).tolist()]
beta_tx_tl = [[0, 0, -5 * tl_scale]]
beta_int_tl = [(-10 * tl_scale * np.array([1.5, 0, 0,
                                           0, 2.0, 0])).tolist()]
sigma_eps_tl = [eps_scale * tl_scale]

################################################################################
# DP sub-model
################################################################################

gamma_tl_dp = [dp_scale]
log_alpha_dp = [np.log(3)]

################################################################################
# SAE sub-model
################################################################################

log_alpha_sae = [np.log(10)]
beta_tx_sae = [(sae_scale * np.array([1, 3, 5])).tolist()]
beta_int_sae = [(sae_scale * np.array([1, 1, 1, 0, 0, 0])).tolist()]

################################################################################
# PS sub-model
################################################################################

beta_0_ps = [1]
beta_1_ps = [0]
beta_tx_ps = [(-10 * ps_scale * np.array([1.5, 1, 1])).tolist()]
beta_tl_ps = [-5e0 * ps_scale / tl_scale]
beta_sae_ps = [(-ps_scale * np.array([0.5, 1, 1.25])).tolist()]
sigma_eps_ps = [10 * eps_scale * ps_scale]

################################################################################
# Patient-level noise
################################################################################

Omega_u = [np.eye(5).tolist()]
sigma_u = [(u_scale * np.array([10 * tl_scale,
                                5 * tl_scale,
                                ps_scale,
                                ps_scale,
                                4 * sae_scale])).tolist()]
      
inputs = dict(
    # data size
    N_pt = N_pt,
    N_bm = N_bm,
    N_tx = N_tx,
    N_ps = N_ps,
    N_tl = N_tl,
    N_samples = N_samples,

    # data structure
    bm_inputs = bm_inputs,
    tx_indicators = tx_indicators,
    x_loc = x_loc,
    t_tl = t_tl,
    pt_idx_tl = pt_idx_tl,
    t_ps = t_ps,
    pt_idx_ps = pt_idx_ps,

    # TL parameters
    beta_0_tl = beta_0_tl,
    beta_1_tl = beta_1_tl,
    beta_bm_tl = beta_bm_tl,
    beta_tx_tl = beta_tx_tl,
    beta_int_tl = beta_int_tl,
    sigma_eps_tl = sigma_eps_tl,

    # DP parameters
    gamma_tl_dp = gamma_tl_dp,
    log_alpha_dp = log_alpha_dp,

    # SAE parameters
    beta_tx_sae = beta_tx_sae,
    beta_int_sae = beta_int_sae,
    log_alpha_sae = log_alpha_sae,

    # PS parameters
    beta_0_ps = beta_0_ps,
    beta_1_ps = beta_1_ps,
    beta_tx_ps = beta_tx_ps,
    beta_tl_ps = beta_tl_ps,
    beta_sae_ps = beta_sae_ps,
    sigma_eps_ps = sigma_eps_ps,

    # patient-level parameters
    Omega_u = Omega_u,
    sigma_u = sigma_u,
)

with open("predict_inputs.json", "w") as f:
    json.dump(inputs, f)
