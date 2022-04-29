"""
Copyright 2019 xCures, Inc.
All rights reserved.
All open source components are property of the respective authors.

Author: Asher Wasserman
Email: awasserman@xcures.com
"""

import numpy as np

def random(y_hat):
    """
    Randomly choose a treatment.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples)-shaped array of predicted outcomes

    Returns 
    -------
    tx : int
       Chosen treatment index
    """
    return np.random.choice(y_hat.shape[0])


def thompson(y_hat):
    """
    Choose a treatment via Thompson sampling.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples)-shaped array of predicted outcomes

    Returns 
    -------
    tx : int
       Chosen treatment index
    """
    y = y_hat[:, np.random.choice(y_hat.shape[1])]
    return np.argmax(y)


def selfish(y_hat):
    """
    Choose a treatment with the highest median.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples)-shaped array of predicted outcomes

    Returns 
    -------
    tx : int
       Chosen treatment index
    """
    return np.argmax(np.median(y_hat, axis=1))


def selfless(y_hat):
    """
    Choose a treatment with the highest variance.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples)-shaped array of predicted outcomes

    Returns 
    -------
    tx : int
       Chosen treatment index
    """
    return np.argmax(np.std(y_hat, axis=1))


def balanced(y_hat):
    """
    If there is an obviously better treatment, choose that.
    If there is a tie, choose the one with the highest variance.
    "not obviously worse" here means the the lower median is greater than the 
    higher median minus the higher one's standard deviation.

    Parameters
    ----------
    y_hat : ndarray
        (N_tx, N_samples)-shaped array of predicted outcomes

    Returns 
    -------
    tx : int
       Chosen treatment index
    """
    N_tx = y_hat.shape[0]
    med = np.median(y_hat, axis=1)
    std = np.std(y_hat, axis=1)

    imax = np.argmax(med)
    best_treatments = [imax]
    best_med = med[imax]
    best_std = std[imax]
    for i in range(N_tx):
        if i == imax:
            continue
        if med[i] > best_med - best_std:
            best_treatments.append(i)
    if len(best_treatments) == 1:
        return imax
    else:
        best_stds = [std[i] for i in best_treatments]
        return best_treatments[np.argmax(best_stds)]
