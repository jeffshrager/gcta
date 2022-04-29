"""
Copyright 2019 xCures, Inc.
All rights reserved.
All open source components are property of the respective authors.

Author: Asher Wasserman
Email: awasserman@xcures.com
"""

import os
from os import path
import pickle
import logging

from ordered_set import OrderedSet
import numpy as np
import pandas as pd
import pystan

from .data import PatientData

standir = path.join(path.abspath(path.dirname(__file__)), "stan")

class suppress_stdout_stderr:
    """
    Suppressing stan output is disappointingly difficult:
    https://github.com/facebook/prophet/issues/223

    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    """
    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for i in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *args):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        for fd in self.null_fds + self.save_fds:
            os.close(fd)
    

def cmdstan_to_pystan_columns(column_name):
    if "[" not in column_name:
        return column_name
    else:
        return column_name.replace("]", "").replace("[", ".")

    
def sample_nuts(sm, inputs, pars=None, chains=1, warmup=200,
                size=1000, verbose=False, **kwargs):
    if verbose:
        fit = sm.sampling(data=inputs, pars=pars, chains=chains, warmup=warmup,
                          iter=size + warmup, **kwargs)
    else:
        with suppress_stdout_stderr():
            fit = sm.sampling(data=inputs, pars=pars, chains=chains, warmup=warmup,
                              iter=size + warmup, **kwargs)
    df = fit.to_dataframe(pars=pars)
    # convert pystan column names to cmdstan column names
    df.columns = list(map(cmdstan_to_pystan_columns, df.columns))
    return df


def sample_advi(sm, inputs, pars=None, size=1000, verbose=False,
                algorithm="meanfield", fail_fast=False, **kwargs):
    try:
        if verbose:
            res = sm.vb(data=inputs, pars=pars, output_samples=size,
                        algorithm=algorithm, **kwargs)
        else:
            with suppress_stdout_stderr():
                res = sm.vb(data=inputs, pars=pars, iter=size)
    except RuntimeError as err:
        if fail_fast:
            raise RuntimeError("ADVI failed: " + err)
        logging.warning("ADVI failed, failing back on NUTS")
        return sample_nuts(sm, inputs, pars=pars, size=size, verbose=verbose,
                           **kwargs)
    outfile = res["args"]["sample_file"].decode()
    df = pd.read_csv(outfile, comment="#")
    os.remove(outfile)
    # convert pystan column names to cmdstan column names
    df.columns = list(map(cmdstan_to_pystan_columns, df.columns))
    # oddly enough, it draws (size + 1) samples...
    df = df.iloc[1:, :]
    return df


class Model:

    def forward(self, **kwargs):
        """Predict data from model parameters"""
        raise NotImplementedError

    def backward(self, **kwargs):
        """Infer model parameters from data"""
        raise NotImplementedError

    
class StanModel(Model):

    def __init__(self, prefix):
        """
        Parameters
        ----------
        prefix : str
            prefix for model code
        """
        self._prefix = prefix
        forward_prefix = prefix + "_predict"
        backward_prefix = prefix + "_infer"
        prior_prefix = prefix + "_prior"
        self.forward_model = self.load_stan_model(forward_prefix)
        self.backward_model = self.load_stan_model(backward_prefix)
        self.prior_model = self.load_stan_model(prior_prefix)

    def load_stan_model(self, prefix):
        pkl_file = path.join(standir, prefix + ".pkl")
        stan_file = path.join(standir, prefix + ".stan")        
        try:
            with open(pkl_file, "rb") as f:
                stan_model = pickle.load(f)
        except FileNotFoundError:
            stan_model = pystan.StanModel(stan_file)
            with open(pkl_file, "wb") as f:
                pickle.dump(model, f)
        return stan_model


class LinearStanModel(StanModel):
    """Linear response Stan Model"""
    def __init__(self, prefix, N_bm, N_tx, hyperparameters, pars, texs, sizes):
        """
        Parameters
        ----------
        N_bm : int
            number of biomarkers
        N_tx : int
            number of treatments
        hyperparameters : dict
            dictionary of hyperparameters
        """
        self.N_bm = N_bm
        self.N_tx = N_tx
        self.N_int = N_bm * N_tx
        self.hyperparameters = hyperparameters
        columns = []
        labels = []
        for j, name in enumerate(pars):
            tex = texs[j]
            n = sizes[j]
            if n == 1 and (name not in ["beta_bm", "beta_tx"]):
                columns.append(name)
                labels.append(tex)
            else:
                for i in range(1, n + 1):
                    columns.append(f"{name}.{i}")
                    labels.append(tex.format(i=i))
        self._columns = columns
        self._labels = labels
        self._pars = pars
        self._par_sizes = sizes
        self.N_columns = sum(sizes)
        super().__init__(prefix=prefix)

    def __repr__(self):
        return f"<{self.__class__.__name__}: N_bm = {self.N_bm}, N_tx = {self.N_tx}>"

    def forward(self, theta, x, tx, t, size=1, warmup=200, verbose=False,
                **kwargs):
        """
        Parameters
        ----------
        theta : ndarray
            (size, N_params)-sized array of parameter samples
        x : array_like
            patient biomarker indicators
        tx : int
            treatment index (0-based)
        t : float or array_like
            time (in months) for predictions
        size : int, optional
            number of samples to draw, default to 1
        warmup : int, optional
            number of warmup iterations, default to 200

        Returns
        -------
        y : ndarray
            (size, len(t))-sized array of predicted outcomes
        """
        inputs = self.construct_forward_inputs(theta, x, tx, t)
        df = sample_nuts(self.forward_model, inputs=inputs, warmup=warmup,
                         size=size, verbose=verbose, **kwargs)
        pars = [f"y_res.{i + 1}" for i in range(inputs["N_res"])]
        return df.loc[:, pars].values
            
    def backward(self, data=None, size=1000, warmup=200, use_vb=False,
                 verbose=False, **kwargs):
        """
        Parameters
        ----------
        data : PatientData instance
        size : int, optional
            number of samples to draw, default to 1000
        warmup : int, optional
            number of warmup iterations, default to 200

        Returns
        -------
        theta : ndarray
            (size, N_params)-sized array of parameter samples
        """
        if data is None:
            data = PatientData(self.N_bm, self.N_tx)
        inputs = self.construct_backward_inputs(data)
        if len(data) == 0:
            # sample from prior
            sm = self.prior_model
        else:
            sm = self.backward_model
        if use_vb:
            df = sample_advi(sm, inputs, size=size, verbose=verbose, **kwargs)
        else:
            df = sample_nuts(sm, inputs, size=size, warmup=warmup,
                             verbose=verbose, **kwargs)
        return df.loc[:, self._columns].values

    def construct_forward_inputs(self, theta, x, tx, t):
        """
        Parameters
        ----------
        theta : ndarray
            (size, N_params)-sized array of parameter samples
        x : array_like
            patient biomarker indicators
        tx : int
            treatment index (0-based)
        t : array_like
            time (in months) for predictions

        Returns
        -------
        inputs : dict
            dictionary of inputs for Stan
        """
        try:
            N_res = len(t)
        except TypeError:
            N_res = 1
            t = [t]
        t_res = np.array(t)            
        N_pt = 1
        N_bm = self.N_bm
        N_tx = self.N_tx
        a = np.zeros(N_tx).astype(int)
        a[tx] = 1
        assert N_bm == len(x)
        N_int = N_bm * N_tx
        bm_indicators = np.array(x).reshape((1, N_bm)).astype(int).tolist()
        tx_indicators = np.array(a).reshape((1, N_tx)).astype(int).tolist()
        if len(theta.shape) == 1:
            theta = theta.reshape((1, theta.size))
        N_samples = theta.shape[0]
        pt_index = [1 for i in range(N_res)]
        inputs = {}
        inputs["N_samples"] = N_samples
        inputs["N_bm"] = N_bm
        inputs["N_tx"] = N_tx
        inputs["N_pt"] = N_pt
        inputs["N_res"] = N_res
        inputs["bm_indicators"] = bm_indicators
        inputs["tx_indicators"] = tx_indicators
        inputs["t_res"] = t_res
        inputs["pt_index"] = pt_index
        start = 0
        for name, size in zip(self._pars, self._par_sizes):
            end = start + size
            if size == 1:
                if name not in ["beta_bm", "beta_tx"]:
                    inputs[name] = theta[:, start]
                else:
                    # handle single biomarker/treatment model
                    inputs[name] = theta[:, start].reshape((-1, 1))
            else:
                inputs[name] = theta[:, start:end]
            start = end    
        return inputs

    def construct_backward_inputs(self, data):
        """
        Parameters
        ----------
        data : PatientData instance

        Returns
        -------
        inputs : dict
            dictionary of inputs for Stan
        """
        biomarker_data = data.biomarkers
        patient_indices = data.biomarkers.index
        df = data.outcomes
        N_bm = self.N_bm
        N_tx = self.N_tx
        # Make N_pt "pseudopatients", one per observed treatment per real patient
        N_pt_orig = len(df.index.levels[0])
        N_tx_for_pt = np.array([len(set(df.loc[pt].index.codes[0]))
                                for pt in patient_indices])
        N_pt = int(np.sum(N_tx_for_pt))
        bm_indicators = np.zeros((N_pt, N_bm))
        tx_indicators = np.zeros((N_pt, N_tx))
        i = 0
        for pt in patient_indices:
            for tx in OrderedSet(df.loc[pt].index.codes[0]):
                bm_indicators[i] = biomarker_data.loc[pt]
                tmp = np.zeros(N_tx)
                tmp[tx] = 1
                tx_indicators[i] = tmp
                i += 1
        # construct response vector
        N_res = df.shape[0]
        t_res = np.array([df.index.levels[2][i] for i in df.index.codes[2]])
        y_res = df.y.values
        # add one to account for stan's 1-based indexing
        pt_index = np.array([df.index.levels[0][i] for i in df.index.codes[0]]) + 1
        inputs = {}
        inputs["N_bm"] = N_bm
        inputs["N_tx"] = N_tx
        if len(data) > 0:
            inputs["N_pt"] = N_pt
            inputs["N_res"] = N_res
            inputs["bm_indicators"] = bm_indicators.astype(int)
            inputs["tx_indicators"] = tx_indicators.astype(int)
            inputs["t_res"] = t_res
            inputs["pt_index"] = pt_index.astype(int)
            inputs["y_res"] = y_res
        inputs.update(self.hyperparameters)
        return inputs
    
    
class LMStanModel(LinearStanModel):
    """Linear Multi-level Stan Model"""
    def __init__(self, N_bm, N_tx, hyperparameters):
        pars = ("beta0", "beta_bm", "beta_tx", "beta_int", "mu_u1", "tau_u",
                "rho_u", "sigma_eps")
        texs = ("$\\beta_0$", "$\\beta_\\mathrm{{bm, {i}}}$",
                "$\\beta_\\mathrm{{tx, {i}}}$",
                "$\\beta_\\mathrm{{int, {i}}}$", "$\\mu_\\mathrm{{u,1}}$",
                "$\\tau_{{u, {i}}}$", "$\\rho_u$", "$\\sigma_\\epsilon$")
        sizes = (1, N_bm, N_tx, N_bm * N_tx, 1, 2, 1, 1)
        super().__init__(prefix="lm", N_bm=N_bm, N_tx=N_tx,
                         hyperparameters=hyperparameters,
                         pars=pars, texs=texs, sizes=sizes)
    
                                 
class LFStanModel(LinearStanModel):
    """Linear Fixed-effects Stan Model"""
    def __init__(self, N_bm, N_tx, hyperparameters):
        pars = ("beta0", "beta_bm", "beta_tx", "beta_int", "sigma_eps")
        texs = ("$\\beta_0$", "$\\beta_\\mathrm{{bm, {i}}}$",
                "$\\beta_\\mathrm{{tx, {i}}}$",
                "$\\beta_\\mathrm{{int, {i}}}$", "$\\sigma_\\epsilon$")
        sizes = (1, N_bm, N_tx, N_bm * N_tx, 1)
        super().__init__(prefix="lf", N_bm=N_bm, N_tx=N_tx,
                         hyperparameters=hyperparameters,
                         pars=pars, texs=texs, sizes=sizes)
    
                                 
