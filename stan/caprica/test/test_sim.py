import json
import numpy as np

from tqdm import tqdm

from caprica import sim, policies, plots
from caprica.models import LMStanModel
from caprica.data import PatientData

#####################
# simulation workflow
#####################

hyperparameters = json.load(open("hp.json"))
N_bm = 5
N_tx = 5
N_pt = 50
N_samples = 1000

bm_weights = np.array([0.2, 0.07, 0.44, 0.15, 0.12])

# true data generating model
theta_true = np.array([
    50, # beta0
    -5, -5, -5, -5, -5, # beta_bm
    -10, -5, -5, 0, -5, # beta_tx
    20, 0, 0, 0, 0,     # beta_int,1
    0, 15, 0, 0, 0,     # beta_int,2
    0, 0, 30, 0, 0,     # beta_int,3
    0, 0, 0, 15, 0,     # beta_int,4
    0, 0, 0, 0, 25,     # beta_int,5
    0, 5, 5, 0.2, 5  # mu_u1, tau_u0, tau_u1, rho_u, sigma_eps
])

data = PatientData(N_bm=N_bm, N_tx=N_tx)
model = LMStanModel(N_bm=N_bm, N_tx=N_tx, hyperparameters=hyperparameters)
policy = policies.random

# sample sets of patient biomarkers
x_pt = np.array([sim.biomarker_selector(bm_weights) for i in range(N_pt)])

# sample from prior
theta_prior = sim.infer(model=model, data=data, size=N_samples)

# sample response at t = 0
y0 = np.array([sim.predict(model=model, theta=theta_true, x=x, treatments=0,
                           t=0, size=1)
               for x in x_pt])
# set time difference
dt = 1

thetas = np.zeros((N_pt + 1, N_samples, theta.shape[1]))
thetas[0] = theta_prior

for i, x in tqdm(enumerate(x_pt), total=N_pt):

    # predict patient response for the pt-th patient
    y_hat = sim.predict(model=model, x=x, theta=thetas[i], t=dt)

    # decide on a treatment
    tx = sim.decide(policy, y_hat)
    
    # append initial patient data
    data.append(inplace=True, pt=i + 1, x=x, y=y0[i], tx=tx, t=0)

    # observe actual patient response
    y_obs = sim.observe(model=model, theta_true=theta_true, x=x, tx=tx)
    
    # append observed patient data after dt
    data.append(inplace=True, pt=i + 1, y=y_obs, tx=tx, t=dt)
    
    # sample posterior
    thetas[i + 1] = sim.infer(model=model, data=data, use_vb=True, size=N_samples)
    
# do actual MCMC sampling with final dataset
thetas[-1] = sim.infer(model=model, data=data, use_vb=False, verbose=True,
                       size=N_samples)
    
    
