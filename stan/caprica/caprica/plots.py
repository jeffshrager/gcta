"""
Copyright 2019 xCures, Inc.
All rights reserved.
All open source components are property of the respective authors.

Author: Asher Wasserman
Email: awasserman@xcures.com
"""
from itertools import cycle

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from corner import corner

from .data import PatientData
from .models import LMStanModel
from . import sim

sns.set(style="ticks")

def parse_params(n_params, model=None, parameters=None):
    if parameters is None:
        indices = np.arange(n_params)
    else:
        if all([isinstance(p, str) for p in parameters]):
            assert model is not None, "need to pass in model for named parameters!"
            indices = [model._columns.index(p) for p in parameters]
        else:
            assert all([isinstance(p, int) for p in parameters])
            indices = parameters 
        n_params = len(parameters)
    if model is None:
        labels = [f"p{i}" for i in range(n_params)]
    else:
        labels = [model._labels[i] for i in indices]
    return indices, labels


def traceplot(theta, n_burn=0, theta_true=None, theta_prior=None, model=None,
              parameters=None):
    """Plot parameter values vs MCMC iteration index and show 1D marginalized
    histograms.  Should only be done for MCMC samples, not VB samples.
    
    Parameters
    ----------
    theta : ndarray
        (N_samples, N_columns)-shaped array of posterior samples
    n_burn : int, optional
        Number of iterations to discard
    theta_true : array_like, optional
        True parameters for comparison
    theta_prior : ndarray, optional
        (N_samples, N_columns)-shaped array of prior samples
    model : Model instance, optional
    parameters : list, optional
        list of parameter indices or names (needs model)

    Returns
    -------
    fig : matplotlib.figure.Figure instance
    """
    n_samples = theta.shape[0] - n_burn
    indices, labels = parse_params(n_params=theta.shape[1], model=model,
                                   parameters=parameters)
    x = theta[:, indices]
    
    n_params = len(labels)
    fig, axes = plt.subplots(ncols=2, nrows=n_params, sharey="row",
                             figsize=(9, 3 * n_params),
                             gridspec_kw=dict(width_ratios=[2, 1], wspace=0.1,
                                              hspace=0.3))

    it = np.arange(n_burn, n_samples + n_burn) + 1
    for i, label in enumerate(labels):
        ax = axes[i, 0]
        y = x[n_burn:, i]
        ax.plot(it, y, color="b")
        ax.set_xlabel("MCMC Iteration")
        ax.set_ylabel(label)

        ax = axes[i, 1]
        ax.hist(y, color="b", density=True, bins="auto",
                histtype="stepfilled", alpha=0.5,
                orientation="horizontal", label="Posterior")
        if theta_prior is not None:
            y_prior = theta_prior[:, indices[i]]
            ax.hist(y_prior, color="g", density=True, bins="auto",
                    histtype="step", lw=2,
                    orientation="horizontal", label="Prior")
        if theta_true is not None:
            y_true = theta_true[indices[i]]
            ax.axhline(y_true, color="r", ls="--", label="Truth")
        else:
            y_true = None
        ax.set_xticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        dy = np.ptp(y)
        ymin = np.amin(y) - 0.1 * dy
        ymax = np.amax(y) + 0.1 * dy
        if y_true is not None:
            ymin = min([ymin, y_true])
            ymax = max([ymax, y_true])
        ax.set_ylim(ymin, ymax)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
    return fig


def cornerplot(theta, n_burn=0, theta_true=None, model=None, parameters=None):
    """Plot 2D marginalized histograms.

    Parameters
    ----------
    theta : ndarray
        (N_samples, N_columns)-shaped array of posterior samples
    n_burn : int, optional
        Number of iterations to discard
    theta_true : array_like, optional
        True parameters for comparison
    model : Model instance, optional
    parameters : list, optional
        list of parameter names as strings

    Returns
    -------
    fig : matplotlib.figure.Figure instance
    """
    indices, labels = parse_params(n_params=theta.shape[1], model=model,
                                   parameters=parameters)
    if theta_true is None:
        truths = None
    else:
        truths = theta_true[indices]
    
    x = theta[n_burn:, indices]
    fig = corner(x, labels=labels, truths=truths, truth_color="r", smooth=0.1,
                 show_titles=True, plot_datapoints=False, plot_density=False,
                 color="b")
    return fig


def simplot(thetas, theta_true=None, theta_prior=None, model=None,
            parameters=None):
    """Plot the change in posteriors over time.

    Parameters
    ----------
    thetas : ndarray
        (N_pt, N_samples, N_columns)-shaped array of posterior samples after 
        each patient
    theta_true : array_like, optional
        True parameters for comparison
    theta_prior : ndarray, optional
        (N_samples, N_columns)-shaped array of prior samples
    model : Model instance, optional
    parameters : list, optional
        list of parameter names as strings

    Returns
    -------
    fig : matplotlib.figure.Figure instance
    """
    indices, labels = parse_params(n_params=thetas.shape[2], model=model,
                                   parameters=parameters)
    
    xlow, xmed, xhigh  = np.percentile(thetas[:, :, indices], q=[16, 50, 84],
                                       axis=1)
    # need to transpose because numpy is silly
    x = thetas[-1, :, indices].T
    
    n_params = len(labels)
    n_pt = thetas.shape[0]
    fig, axes = plt.subplots(ncols=2, nrows=n_params, sharey="row",
                             figsize=(9, 3 * n_params),
                             gridspec_kw=dict(width_ratios=[2, 1], wspace=0.1,
                                              hspace=0.3))
    if model is None:
        N_bm = 5
    else:
        N_bm = model.N_bm
    t = np.arange(0, n_pt)
    for i, label in enumerate(labels):
        ax = axes[i, 0]
        ymed = xmed[:, i]
        ylow = xlow[:, i]
        yhigh = xhigh[:, i]        
        ax.plot(t, ymed, color="b")
        ax.fill_between(t, ylow, yhigh, color="b", alpha=0.2)
        ax.set_xlabel("Number of patients")
        ax.set_ylabel(label)

        ax = axes[i, 1]
        y = x[:, i]
        ax.hist(y, color="b", density=True, bins="auto",
                histtype="stepfilled", alpha=0.5,
                orientation="horizontal", label="Posterior")
        if theta_prior is not None:
            y_prior = theta_prior[:, indices[i]]
            ax.hist(y_prior, color="g", density=True, bins="auto",
                    histtype="step", lw=2,
                    orientation="horizontal", label="Prior")
        if theta_true is not None:
            y_true = theta_true[indices[i]]
            ax.axhline(y_true, color="r", ls="--", label="Truth")
        else:
            y_true = None
        ax.set_xticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        slc = slice(2 * N_bm, None)
        dy = np.amax(yhigh[slc] - ylow[slc])
        ymin = np.amin(ylow[slc]) - 0.1 * dy
        ymax = np.amax(yhigh[slc]) + 0.1 * dy
        if y_true is not None:
            ymin = min([ymin, y_true])
            ymax = max([ymax, y_true])
        ax.set_ylim(ymin, ymax)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
    return fig

    
def ppcplot(y_hat, times=None, data=None):
    """Posterior predictive checks.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples, N_times)-shaped array of predicted patient outcomes
    times : array_like, optional
        times (months after treatment) of predictions, if None, defaults to 1
    data : PatientData instance, optional
        if None, just plot the predictions
    """
    N_tx = y_hat.shape[0]
    fig, axes = plt.subplots(nrows=N_tx,
                             figsize=(6, 3 * N_tx),
                             gridspec_kw=dict(hspace=0.3))
    if N_tx == 1:
        axes = [axes]
    if times is None:
        times = np.arange(y_hat.shape[2])
    color = cycle([f"C{i}" for i in range(5)])
    for tx, ax in enumerate(axes):
        ylow, ymed, yhigh = np.percentile(y_hat[tx], axis=0, q=[16, 50, 84])
        c = next(color)
        ax.plot(times, ymed, color=c, ls="-")
        ax.fill_between(times, ylow, yhigh, color=c, alpha=0.2)
        ax.set_ylabel(f"Treatment {tx} response")
        ax.set_xlabel("Time (months)")
    return fig

