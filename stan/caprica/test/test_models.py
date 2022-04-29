import json
import pickle
import numpy as np

from caprica.models import LMStanModel
from caprica.data import PatientData

hyperparameters = json.load(open("hp.json"))

N_tx = 5
N_bm = 5

with open("data.pkl", "rb") as f:
    data = pickle.load(f)
    
model = LMStanModel(N_tx, N_bm, hyperparameters)

x = [0, 1, 1, 0, 0]
tx = 2
t = [1]
theta = model.backward(data, verbose=True, warmup=200, size=500)
y = model.forward(np.median(theta, axis=0), x=x, tx=tx, t=t, verbose=True, size=1)

