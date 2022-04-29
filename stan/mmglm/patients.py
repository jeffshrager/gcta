import numpy as np

randn = np.random.randn

tmin = 0
tmax = 10
dt = 0.5
t = np.arange(tmin, tmax, dt)
N_data_per_pt = t.size

N_bm = 2
N_tx = 3
N_rep = 30
N_pt = N_rep * N_bm**2 * N_tx
N_int = N_bm * N_tx

N_ps = N_pt * N_data_per_pt
N_tl = N_pt * N_data_per_pt

x_unique = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
bm_inputs = np.tile(x_unique, (N_rep * N_tx, 1)).tolist()

tx_unique = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
tx_indicators = np.zeros((N_pt, N_tx))

tx_indicators = np.tile(np.array([[1, 0, 0],
                                  [1, 0, 0],
                                  [1, 0, 0],
                                  [1, 0, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [0, 1, 0],
                                  [0, 0, 1],
                                  [0, 0, 1],
                                  [0, 0, 1],
                                  [0, 0, 1]]), (N_rep, 1)).tolist()

x_loc = [1 for i in range(N_pt)]

t_tl = np.vstack([t for i in range(N_pt)]).flatten().tolist()
pt_idx_tl = np.vstack([i * np.ones(N_data_per_pt)
                       for i in range(1, 1 + N_pt)]).flatten().astype(int).tolist()
t_ps = np.vstack([t for i in range(N_pt)]).flatten().tolist()
pt_idx_ps = np.vstack([i * np.ones(N_data_per_pt)
                        for i in range(1, 1 + N_pt)]).flatten().astype(int).tolist()
