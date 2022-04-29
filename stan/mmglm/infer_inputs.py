import json
import numpy as np
import pandas as pd

with open("patients.py") as f:
    exec(f.read())

# load data
df = pd.read_csv("predict_outputs.csv", index_col=0)
y_tl = [df.loc[0, f"y_tl[{pt}]"] for pt in range(1, N_tl + 1)]
y_ps = [df.loc[0, f"y_ps[{pt}]"] for pt in range(1, N_ps + 1)]
T_dp = [df.loc[0, f"T_dp[{pt}]"] for pt in range(1, N_pt + 1)]
T_sae = [df.loc[0, f"T_sae[{pt}]"] for pt in range(1, N_pt + 1)]

# assume data are uncensored
censored_dp = [0 for i in range(N_pt)]
censored_sae = [0 for i in range(N_pt)]

time_scale = 6    # months
tl_scale = 1e-2   # delta TL per month
dp_scale = 1e-3   # delta DP hazard per TL per month
sae_scale = 1e-7  # delta SAE hazard per month
ps_scale = 1e-2   # delta PS per month
eps_scale = 1     # S/N ratio
u_scale = 1       # S/N ratio

data = dict(
    # data size
    N_pt = N_pt,
    N_bm = N_bm,
    N_tx = N_tx,
    N_ps = N_ps,
    N_tl = N_tl,

    # data structure
    bm_inputs = bm_inputs,
    tx_indicators = tx_indicators,
    x_loc = x_loc,
    t_tl = t_tl,
    pt_idx_tl = pt_idx_tl,
    t_ps = t_ps,
    pt_idx_ps = pt_idx_ps,

    # response data
    y_tl = y_tl,
    y_ps = y_ps,
    T_dp = T_dp,
    censored_dp = censored_dp,
    T_sae = T_sae,
    censored_sae = censored_sae,
    
    # TL hyperparameters
    loc_0_tl = 0,
    scale_0_tl = time_scale * tl_scale,
    loc_1_tl = 0,
    scale_1_tl = tl_scale,
    loc_bm_tl = np.zeros(N_bm).tolist(),
    scale_bm_tl = (tl_scale * np.ones(N_bm)).tolist(),
    loc_tx_tl = np.zeros(N_tx).tolist(),
    scale_tx_tl = (tl_scale * np.ones(N_tx)).tolist(),
    loc_int_tl = np.zeros(N_int).tolist(),
    scale_int_tl = (tl_scale * np.ones(N_int)).tolist(),
    scale_eps_tl = eps_scale * tl_scale,
    
    # DP hyperparameters
    loc_gamma_tl_dp = 0,
    scale_gamma_tl_dp = dp_scale,
    loc_log_alpha_dp = 0,
    scale_log_alpha_dp = 0.1,

    # SAE hyperparameters
    loc_tx_sae = np.zeros(N_tx).tolist(),
    scale_tx_sae = (sae_scale * np.ones(N_tx)).tolist(),
    loc_int_sae = np.zeros(N_int).tolist(),
    scale_int_sae = (sae_scale * np.ones(N_int)).tolist(),
    loc_log_alpha_sae = 0,
    scale_log_alpha_sae = 0.1,

    # PS hyperparameters
    loc_0_ps = 0,
    scale_0_ps = time_scale * ps_scale,
    loc_1_ps = 0,
    scale_1_ps = ps_scale,
    loc_tx_ps = np.zeros(N_tx).tolist(),
    scale_tx_ps = (ps_scale * np.ones(N_tx)).tolist(),
    loc_sae_ps = np.zeros(N_tx).tolist(),
    scale_sae_ps = (ps_scale * np.ones(N_tx)).tolist(),
    scale_eps_ps = eps_scale * ps_scale,

    # Patient-level hyperparameters
    eta_u = 2,
    # scale_u = (u_scale * np.array([tl_scale, tl_scale, sae_scale, ps_scale,
    #                                ps_scale])).tolist()
    scale_u = (u_scale * np.array([tl_scale, tl_scale])).tolist()
)

json.dump(data, open("infer_inputs.json", "w"))
