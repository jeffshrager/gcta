from glob import glob
import re
from collections import namedtuple
from itertools import cycle
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from corner import corner

Parameter = namedtuple("Parameter", 
                       ["name", "tex", "prior", "value"])

def normal(x, mean, std):
    return stats.norm.pdf(x, loc=mean, scale=std)

def hcauchy(x, scale):
    return stats.cauchy.pdf(x, loc=0, scale=scale)

def beta(x, a):
    return stats.beta.pdf(x, a, a)


def get_truths(predict_input):
    """Extract a dictionary of true parameter values from the predict json input.
    """
    with open(predict_input) as f:
        data = json.loads(f.read())
    names = ("beta0", "beta_bm", "beta_tx", "beta_int", 
             "tau_u", "mu_u1", "rho_u", "sigma_eps")
    truth_dict = {}
    for name in names:
        truth = data[name][0]
        try:
            size = len(truth)
        except TypeError:
            truth_dict[name] = truth
            continue
        for k, entry in enumerate(truth):
            truth_dict[name + f".{k + 1}"] = entry
    return truth_dict


def make_parameters(infer_input=None, predict_input=None, N_bm=5, N_tx=5):

    # make names and priors
    if infer_input is not None:
        with open(infer_input) as f:
            data = json.loads(f.read())
        N_bm = data["N_bm"]
        N_tx = data["N_tx"]
        N_int = N_bm * N_tx
        priors = {
            "mu_u1": lambda x: normal(x, data["loc_mu_u1"], data["scale_mu_u1"]),
            "tau_u.1": lambda x: hcauchy(x, data["scale_tau_u"][0]),
            "tau_u.2": lambda x: hcauchy(x, data["scale_tau_u"][1]),
            "rho_u": lambda x: beta((x + 1) / 2, data["eta_u"]),
            "sigma_eps": lambda x: hcauchy(x, data["scale_sigma_eps"]),
            "beta0": lambda x: np.ones_like(x) / 100,
        }
        priors.update({f"beta_bm.{i + 1}": (lambda x: normal(x, data["loc_bm"][i], data["scale_bm"][i][i])) for i in range(N_bm)})
        priors.update({f"beta_tx.{i + 1}": (lambda x: normal(x, data["loc_tx"][i], data["scale_tx"][i][i])) for i in range(N_tx)})
        priors.update({f"beta_int.{i + 1}": (lambda x: normal(x, data["loc_int"][i], data["scale_int"][i][i])) for i in range(N_int)})
        names = list(priors.keys())
    else:
        try:
            N_int = N_bm * N_tx
        except:
            raise ValueError("Must pass in parameter counts")
        names = ["mu_u1", "tau_u.1", "tau_u.2", "beta0", "rho_u", "sigma_eps"]
        names.extend([f"beta_bm.{i}" for i in range(1, N_bm + 1)])
        names.extend([f"beta_tx.{i}" for i in range(1, N_tx + 1)])
        names.extend([f"beta_int.{i}" for i in range(1, N_int + 1)])        
        priors = {name: None for name in names}

    # make pretty tex strings
    tex_dict = {"mu_u1": r"$\mu_{u_1}$",
                "tau_u.1": r"$\tau_{u_0}$",
                "tau_u.2": r"$\tau_{u_1}$",
                "rho_u": r"$\rho_u$",
                "sigma_eps": r"$\sigma_\epsilon$",
                "beta0": r"$\beta_0$"}
    tex_dict.update({f"beta_bm.{i}": f"$\\beta_\\mathrm{{bm,{i}}}$"
                     for i in range(1, N_bm + 1)})
    tex_dict.update({f"beta_tx.{i}": f"$\\beta_\\mathrm{{tx,{i}}}$"
                     for i in range(1, N_tx + 1)})
    tex_dict.update({f"beta_int.{i}": f"$\\beta_\\mathrm{{int,{i}}}$"
                     for i in range(1, N_int + 1)})

    # get true values
    if predict_input is not None:
        with open(predict_input) as f:
            truth_dict = get_truths(predict_input)
    else:
        truth_dict = {name: None for name in names}

    parameters = [Parameter(name, tex_dict[name], priors[name], truth_dict[name]) for name in names]
    return parameters


def parse_outputs(outputs):
    """Return a list of data frames containing Stan outputs.
    
    Outputs can be any of:
    1. a list of data frame (will return the list)
    2. a single data frame (will return as a single-entry list)
    3. a string with the file name (will open as a df)
    4. a string with a wildcard expressions representing multiple output files
    5. a list of strings containing file names or wildcard expressions
    
    """
    if isinstance(outputs, pd.DataFrame):
        return [outputs]
    elif isinstance(outputs, str):
        if re.match("[\*\?\]\[]", outputs) is None:
            outputs = [pd.read_csv(outputs, comment="#")]
        else:
            outputs = [pd.read_csv(f, comment="#") for f in glob(outputs)]
        return outputs
    else:
        # is already a list
        if isinstance(outputs[0], pd.DataFrame):
            return outputs
        else:
            new_outputs = []
            for i, entry in enumerate(outputs):
                if re.match("[\*\?\]\[]", outputs) is None:
                    new_outputs.append(entry)
                else:
                    new_outputs.extend(glob(entry))
            return [pd.read_csv(f, comment="#") for f in new_outputs]
            

def parse_parameters(parameters, infer_input=None, predict_input=None):
    if isinstance(parameters, str):
        parameters = [parameters]
    if hasattr(parameters[0], "name"):
        return parameters
    param_list = make_parameters(infer_input=infer_input,
                                 predict_input=predict_input)
    params = []
    for p in param_list:
        if p.name in parameters:
            params.append(p)
    return params


def traceplot(outputs, parameters, n_burn=0,
              infer_input=None, predict_input=None):
    
    outputs = parse_outputs(outputs)
    parameters = parse_parameters(parameters,
                                  infer_input=infer_input,
                                  predict_input=predict_input)
    n_chains = len(outputs)
    n_samples = outputs[0].shape[0] - n_burn
    n_params = len(parameters)
    colors = [f"C{i}" for i in range(6)]
 
    fig, axes = plt.subplots(ncols=2, nrows=n_params, sharey="row", figsize=(9, 3 * n_params),
                             gridspec_kw=dict(width_ratios=[2, 1], wspace=0.1))
   
    for i, p in enumerate(parameters):
        
        x = np.array([df.loc[n_burn:, p.name] for df in outputs]).reshape((n_chains, n_samples))
        
        ax = axes[i, 0]
        for j in range(n_chains):
            ax.plot(np.arange(n_burn, n_samples + n_burn) + 1, x[j, :], color=colors[j])
        ax.set_xlabel("Iteration")
        ax.set_ylabel(p.tex)
        
        ax = axes[i, 1]
        for i in range(n_chains):
            ax.hist(x[j, :], color=colors[j],
                    density=True, bins="auto", histtype="step", orientation="horizontal")
        dx = np.ptp(x)
        xmin = np.amin(x) - 0.1 * dx
        xmax = np.amax(x) + 0.1 * dx
        xx = np.linspace(xmin, xmax)
        if p.prior is not None:
            ax.plot(p.prior(xx), xx, color="r", ls="-", lw=2)
        if p.value is not None:
            ax.axhline(p.value, color="r", ls="--")
        ax.set_xticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_ylim(xmin, xmax)

    fig.show()
    return fig
    
    
def cornerplot(outputs, parameters, n_burn=0,
               infer_input=None, predict_input=None):
    df = parse_outputs(outputs)[0]
    parameters = parse_parameters(parameters,
                                  infer_input=infer_input,
                                  predict_input=predict_input)
    samples = df.loc[n_burn:, [p.name for p in parameters]].values
    labels = [p.tex for p in parameters]
    if parameters[0].value is not None:
        truths = [p.value for p in parameters]
    else:
        truths = None
    fig = corner(samples, truths=truths, labels=labels, color="C0",
                 show_titles=True, truth_color="r",
                 plot_datapoints=False, plot_density=False,
                 fill_contours=True)
    return fig


def ppcplot(outputs, infer_input, n_burn=0, patients=5):
    
    df = parse_outputs(outputs)[0]
    with open(infer_input) as f:
        data = json.loads(f.read())
    N_res = data["N_res"]
    N_pt = data["N_pt"]
    y_res = np.array(data["y_res"])
    t_res = np.array(data["t_res"])
    pt_index = np.array(data["pt_index"])
    y_ppc = df.loc[n_burn:, "y_ppc.1":f"y_ppc.{N_res}"].values
    
    if isinstance(patients, int):
        patients = np.random.choice(np.arange(1, N_pt + 1), size=patients)

    color = cycle([f"C{i}" for i in range(5)])
    for pt in patients:
        pick = (pt == pt_index)
        t = t_res[pick]
        y_data = y_res[pick]
        y_model = y_ppc[:, pick]
        low, med, high = np.percentile(y_model, axis=0, q=[16, 50, 84])
        c = next(color)
        plt.plot(t, y_data, c + "o", mec="k", mew=1, zorder=3)
        plt.plot(t, med, c + "-")
        plt.fill_between(t, low, high, color=c, alpha=0.2)
    plt.ylabel("Response")
    plt.xlabel("Time (months)")
    return plt.gcf()


def residualplot(outputs, infer_input, n_burn=0):
    
    df = parse_outputs(outputs)[0]
    with open(infer_input) as f:
        data = json.loads(f.read())
    N_res = data["N_res"]
    N_pt = data["N_pt"]
    y_res = np.array(data["y_res"])
    t_res = np.array(data["t_res"])
    pt_index = np.array(data["pt_index"])
    y_ppc = df.loc[n_burn:, "y_ppc.1":f"y_ppc.{N_res}"].values
    sigma_eps = np.median(df.loc[n_burn:, "sigma_eps"].values)
    
    low, med, high = np.percentile(y_ppc, axis=0, q=[16, 50, 84])
    residuals = med - y_res
    err = (high - low) / 2
    
    fig, axes = plt.subplots(ncols=2, nrows=1, sharey="row", figsize=(9, 3),
                             gridspec_kw=dict(width_ratios=[2, 1], wspace=0.1))
    
    ax = axes[0]
    color = cycle([f"C{i}" for i in range(5)])
    for pt in range(1, N_pt + 1):
        pick = (pt == pt_index)
        t = t_res[pick] + np.random.randn(np.count_nonzero(pick)) / 10
        dy = residuals[pick]
        dy_err = err[pick]
        c = next(color)
        ax.errorbar(t, dy, yerr=dy_err, marker="o", mfc="none", mec=c, ls="-")
    ax.set_ylabel("(Model - Data) Response")
    ax.set_xlabel("Time (months)")
    
    ax = axes[1]
    ax.hist(residuals, orientation="horizontal", histtype="stepfilled",
            density=True, bins="auto", color="k", alpha=0.5,
            label="Residuals")
    yy = np.linspace(residuals.min(), residuals.max(), 100)
    ax.plot(normal(yy, 0, sigma_eps), yy, color="r",
            label=r"$\mathcal{N}(0, \sigma_\epsilon)$")
    ax.legend(loc="upper right", frameon=False)
    ax.set_xticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    fig.tight_layout()
    return fig


    
def forecastplot(outputs, forecast_input, bm_idx=0):
    
    df = parse_outputs(outputs)[0]
    with open(forecast_input) as f:
        data = json.loads(f.read())
        
    N_res = data["N_res"]
    N_pt = data["N_pt"]
    N_tx = data["N_tx"]
    t_res = np.array(data["t_res"])
    pt_index = np.array(data["pt_index"])
    y_res = df.loc[:, "y_res.1":f"y_res.{N_res}"].values
    ylow, ymed, yhigh = np.percentile(y_res, axis=0, q=[16, 50, 84])
    
    color = cycle([f"C{i}" for i in range(N_tx)])

    for pt in range(1 + bm_idx * N_tx, 1 + (bm_idx + 1) * N_tx):
        
        c = next(color)

        pick = (pt == pt_index)
        t = t_res[pick]
        tx = ((pt - 1)% N_tx) + 1
        plt.plot(t, ymed[pick], color=c, ls="-", label=f"Tx = {tx}")
        plt.fill_between(t, ylow[pick], yhigh[pick], color=c, alpha=0.2)
        
    plt.ylabel("Forecasted response")
    plt.xlabel("Time (months)")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    return plt.gcf()
    
