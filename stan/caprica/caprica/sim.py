"""
Copyright 2019 xCures, Inc.
All rights reserved.
All open source components are property of the respective authors.

Author: Asher Wasserman
Email: awasserman@xcures.com
"""

import numpy as np


def biomarker_selector(weights):
    """
    Parameters
    ----------
    weights : array
        N_bm-sized array of probabilities of having biomarker
    """
    x = np.zeros(weights.shape)
    while not np.any(x):
        x = np.array([np.random.choice([0, 1], p=[1 - w, w])
                      for w in weights])
    return x


def predict(model, theta, x, treatments=None, t=1, size=1000, warmup=1000,
            verbose=False, **kwargs):
    """Predict patient outcomes under each treatment.

    Parameters
    ----------
    model : Model instance
    theta : ndarray
        posterior probability samples, shape (N_samples, N_parameters)
    x : array_like
        patient biomarker indicators, size N_bm
    t : float or array_like, optional
        time (in months) after treatment for outcome predictions
        defaults to 1 month
    treatments : list, optional
        list of treatment indices (0-based), if None, defaults to all
    size : int, optional
        number of samples to draw
        defaults to 1000

    Returns
    -------
    y : ndarray
        (N_tx, N_samples, len(t))-shaped array of patient outcomes
        Squeezes last axis out if len(t) == 1
    """
    if treatments is None:
        treatments = np.arange(model.N_tx)
    try:
        N_treatments = len(treatments)
    except TypeError:
        N_treatments = 1
        treatments = [treatments]
    y_hat = np.array([model.forward(x=x, theta=theta, tx=tx, t=t, size=size,
                                    warmup=warmup, verbose=verbose, **kwargs)
                     for tx in treatments])
    return np.squeeze(y_hat)


def decide(policy, y):
    """Choose a treatment based on predicted outcomes.

    Parameters
    ----------
    policy : Policy instance
    y : ndarray
        (N_tx, N_samples)-shaped array of predicted patient outcomes
    """
    return policy(y)


def observe(model, theta_true, x, tx, t=1, warmup=200, **kwargs):
    """Generate an outcome from the true data-generating model.

    Parameters
    ----------
    model : Model instance
    theta_true : array_like
        true parameter values for the simulation
    x : list
        biomarker indicators
    tx : int
        treatment index, 0-based
    t : float or array_like, optional
        time (in months) after treatment for outcome predictions
        defaults to 1 month

    Returns
    -------
    y_obs : float or ndarray
        (len(t),)-shaped array of observed patient outcomes
        Squeezes last axis out if len(t) == 1
    """
    y_obs = model.forward(theta=theta_true, x=x, tx=tx, t=t, size=1,
                          warmup=warmup, **kwargs)
    return np.squeeze(y_obs)


def infer(model, data, use_vb=False, verbose=False, size=1000, warmup=1000,
          **kwargs):
    """Infer the posterior probability distribution from stored data.

    Parameters
    ----------
    model : Model instance
    data : PatientData instance
    size : int, optional
        number of samples to draw

    Returns
    -------
    theta : ndarray
        (size, N_params)-shaped array of samples from posterior 
        probability distribution
    """
    return model.backward(data, use_vb=use_vb, size=size, warmup=warmup,
                          verbose=verbose, **kwargs)

