import pickle
import numpy as np

from caprica.data import PatientData

N_bm = 5
N_tx = 5

data = PatientData(N_bm, N_tx)

for bm in range(N_bm):
    x = np.zeros(N_bm)
    x[bm] = 1
    data.append(pt=bm, x=x.astype(int), inplace=True)

for pt in data.biomarkers.index:
    for tx in range(N_tx):
        data.append(pt=pt, tx=tx, t=0.0, y=50 + 10 * np.random.randn(), inplace=True)
        data.append(pt=pt, tx=tx, t=1.0, y=50 + 10 * np.random.randn(), inplace=True)

with open("data.pkl", "wb") as f:
    pickle.dump(data, f)
        
