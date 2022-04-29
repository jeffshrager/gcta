#!/bin/env python3

import sys
import argparse
import json
import numpy as np
import pandas as pd


def parse():
    parser = argparse.ArgumentParser(description="pipe stan inputs and outputs")
    parser.add_argument("out_json", action="store", help="filename for outputted json")
    subparsers = parser.add_subparsers(help="sub-command help")

    parser_i = subparsers.add_parser("infer",
                                     help="make inference inputs from predicted outputs")
    parser_i.set_defaults(command="infer")
    parser_i.add_argument("in_csv", help="prediction csv output file")
    parser_i.add_argument("in_json", help="prediction json input file")
    parser_i.add_argument("--hp", "--hyperparameters", default=None,
                          help="hyperparameters dictionary from json file")
    
    parser_p = subparsers.add_parser("predict",
                                     help="make prediction inputs from inferred outputs")
    parser_p.set_defaults(command="predict")
    parser_p.add_argument("in_csv", help="inference csv output file")
    parser_p.add_argument("in_json", help="inference json input file")
    parser_p.add_argument("--bm", "--biomarkers", default=None,
                          help="biomarker indicator from json file")
    parser_p.add_argument("--duration", default=6,
                          help="time duration, in months")
    parser_p.add_argument("--period", default=1,
                          help="time between samples, in months")
    return parser


def make_hyperparameters(N_bm, N_tx,
                         loc_mu_u1=0,
                         scale_mu_u1=20,
                         scale_tau_u=(20, 20),
                         eta_u=1,
                         loc_bm=0,
                         scale_bm=20,
                         loc_tx=0,
                         scale_tx=20,
                         loc_int=0,
                         scale_int=20,
                         scale_sigma_eps=20):
    N_int = N_bm * N_tx
    if isinstance(loc_bm, (int, float)):
        loc_bm = (np.ones(N_bm) * loc_bm).tolist()
    if isinstance(scale_bm, (int, float)):
        scale_bm = (np.eye(N_bm) * scale_bm).tolist()
    if isinstance(loc_tx, (int, float)):
        loc_tx = (np.ones(N_tx) * loc_tx).tolist()
    if isinstance(scale_tx, (int, float)):
        scale_tx = (np.eye(N_tx) * scale_tx).tolist()
    if isinstance(loc_int, (int, float)):
        loc_int = (np.ones(N_int) * loc_int).tolist()
    if isinstance(scale_int, (int, float)):
        scale_int = (np.eye(N_int) * scale_int).tolist()
        
    return dict(eta_u=eta_u,
                loc_mu_u1=loc_mu_u1,
                scale_mu_u1=scale_mu_u1,
                scale_tau_u=scale_tau_u,
                loc_bm=loc_bm,
                scale_bm=scale_bm,
                loc_tx=loc_tx,
                scale_tx=scale_tx,
                loc_int=loc_int,
                scale_int=scale_int,
                scale_sigma_eps=scale_sigma_eps)


def infer_from_predict(predict_csv, predict_json,
                       outfile=None,
                       hp_kwargs=None):
    """
    Construct inference input from prediction output for a recovery test.

    predict_csv : str, csv filename
    predict_json : str, json filename
    outfile : str, output filename
    hp_kwargs : dict, hyperparameter keyword arguments
    """

    # read input data for predict
    with open(predict_json) as f:
        data = json.loads(f.read())
    N_res = data["N_res"]
    N_bm = data["N_bm"]
    N_tx = data["N_tx"]
    
    # extract predicted response
    df = pd.read_csv(predict_csv, comment="#")
    y_res = np.array([df.loc[0, f"y_res.{i}"] for i in range(1, N_res + 1)])

    # remove unknown inputs
    for key in ["N_samples", "beta0", "beta_bm", "beta_tx", "beta_int",
                "mu_u1", "tau_u", "rho_u", "sigma_eps"]:
        data.pop(key)
        
    # add mock responses
    data["y_res"] = y_res.tolist()

    # add hyperparameters
    if hp_kwargs is None:
        hp_kwargs = {}
        
    data.update(make_hyperparameters(N_bm, N_tx, **hp_kwargs))
        
    if outfile is not None:
        with open(outfile, "w") as f:
            f.write(json.dumps(data))
    return data


def predict_from_infer(infer_csv, infer_json, outfile=None,
                       bm_indicators=None, time_samples=None, N_samples=None):
    """
    Construct prediction input from inference output for a forecast.
    
    infer_csv : str, csv filename
    infer_json : str, json filename
    outfile : str, output filename
    bm_indicators : array-like, indicator array of patient biomarkers
    time_samples : array-like, grid of times to predict response
    N_samples : int, number of samples to use from inferred parameters
    """
    # read input data for infer
    with open(infer_json) as f:
        data = json.loads(f.read())

    # remove extraneous inputs
    for key in ["y_res", "eta_u", "loc_mu_u1", "scale_mu_u1", "scale_tau_u",
                "loc_bm", "scale_bm", "loc_tx", "scale_tx",
                "loc_int", "scale_int", "scale_sigma_eps"]:
        data.pop(key)

    N_bm = data["N_bm"]
    N_tx = data["N_tx"]
    N_int = N_bm * N_tx
    
    # patient-biomarker-treatment structure
    if time_samples is None:
        # default to monthly cadence for half a year
        time_samples = np.arange(0, 7, 1)
    if bm_indicators is None:
        # default to one patient per biomarker-treatment combination
        N_pt_orig = N_bm        
        N_pt = N_pt_orig * N_tx
        bm_indicators = np.zeros((N_pt, N_bm))
        for bm in range(N_bm):
            bm_indicators[bm * N_tx: (bm + 1) * N_tx, bm] = np.ones(N_tx)
        bm_indicators = bm_indicators.astype(int).tolist()
    else:
        # make a copy of the bm indicators for each treatment
        N_pt_orig = len(bm_indicators)
        N_pt = N_tx * N_pt_orig
        new_bm_indicators = np.zeros((N_pt, N_bm))
        for pt, indicators in enumerate(bm_indicators):
            new_bm_indicators[pt * N_tx: (pt + 1) * N_tx] = np.tile(indicators, (N_tx, 1))
        bm_indicators = new_bm_indicators.astype(int).tolist()
        
    data["N_pt"] = N_pt
    data["N_res"] = N_pt * time_samples.size
    data["t_res"] = np.tile(time_samples, (N_pt,)).tolist()
    data["n_res"] = [time_samples.size for i in range(N_pt)]
    data["pt_index"] = np.hstack([(i + 1) * np.ones(time_samples.size)
                                  for i in range(N_pt)]).astype(int).tolist()
    data["bm_indicators"] = bm_indicators
    # one treatment per patient
    data["tx_indicators"] = np.tile(np.eye(N_tx), (1, N_pt_orig)).T.astype(int).tolist()
        
    # extract inferred parameters
    df = pd.read_csv(infer_csv, comment="#")

    if N_samples is None:
        N_samples = df.shape[0]
        samples = np.arange(N_samples)
    else:
        samples = np.arange(0, df.shape[0], int(df.shape[0] / (N_samples - 1)))
    data["N_samples"] = N_samples
    
    # add scalar parameters
    for key in ["beta0", "mu_u1", "rho_u", "sigma_eps"]:
        data[key] = df.loc[samples, key].values.tolist()

    # add other random-effects parameters
    for key in ["tau_u"]:
        key1 = key + ".1"
        key2 = key + ".2"
        data[key] = list(zip(df.loc[samples, key1].values,
                             df.loc[samples, key2].values))

    # add fixed-effect parameters
    for key, count in zip(["beta_bm", "beta_tx", "beta_int"],
                          [N_bm, N_tx, N_int]):
        data[key] = np.vstack([df.loc[samples, f"{key}.{i}"]
                               for i in range(1, count + 1)]).T.tolist()

    if outfile is not None:
        with open(outfile, "w") as f:
            f.write(json.dumps(data))
    
    return data


def main():
    parser = parse()
    ns = parser.parse_args()
    if len(sys.argv) < 2 or (not hasattr(ns, "command")):
        parser.print_help()
        return
    if ns.command == "infer":
        if ns.hp is not None:
            with open(ns.hp) as f:
                hyperparameters = json.loads(f.read())
        else:
            hyperparameters={}
        data = infer_from_predict(ns.in_csv, ns.in_json,
                                  outfile=ns.out_json,
                                  hp_kwargs=hyperparameters)
    elif ns.command == "predict":
        if ns.bm is not None:
            with open(ns.bm) as f:
                bm_indicators = json.loads(f.read())
        else:
            bm_indicators = None
        time_samples = np.arange(0, ns.duration + ns.period, ns.period)
        data = predict_from_infer(ns.in_csv, ns.in_json,
                                  outfile=ns.out_json,
                                  bm_indicators=bm_indicators,
                                  time_samples=time_samples)
    return data


if __name__ == "__main__":
    main()
