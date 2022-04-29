"""
Copyright 2019 xCures, Inc.
All rights reserved.
All open source components are property of the respective authors.

Author: Asher Wasserman
Email: awasserman@xcures.com
"""

from warnings import warn
import numpy as np
import pandas as pd

class PatientData:
    
    def __init__(self, N_bm, N_tx, outcomes=None, biomarkers=None):
        self.N_bm = N_bm
        self.N_tx = N_tx
        if outcomes is None:
            outcomes = pd.DataFrame(index=pd.MultiIndex(names=["pt", "tx", "t"],
                                                        levels=((), (), ()),
                                                        codes=((), (), ())),
                                    columns=["y"])
                                                                
        if biomarkers is None:
            biomarkers = pd.DataFrame(columns=[f"x{bm}" for bm in range(1, 1+ N_bm)])
        self.outcomes = outcomes
        self.biomarkers = biomarkers
    
    def __repr__(self):
        s = f"<PatientData: {self.N_pt} patients, {self.N_out} responses>"
        return s

    def __len__(self):
        return self.N_out
    
    @property
    def N_pt(self):
        return self.biomarkers.shape[0]

    @property
    def N_out(self):
        return self.outcomes.shape[0]
    
    def append(self, pt=None, tx=None, t=None, y=None, x=None,
               inplace=False, overwrite=False):
        biomarkers = self._append_biomarkers(pt=pt, x=x,
                                             inplace=inplace,
                                             overwrite=overwrite)
        if pt is None:
            pt = biomarkers.index[-1]
        outcomes = self._append_outcomes(pt=pt, tx=tx, t=t, y=y,
                                         inplace=inplace,
                                         overwrite=overwrite)
        if not inplace:
            return PatientData(self.N_bm, self.N_tx, outcomes, biomarkers)

    def _append_biomarkers(self, pt, x, inplace=False, overwrite=False):
        if inplace:
            biomarkers = self.biomarkers
        else:
            biomarkers = self.biomarkers.copy()
        if pt is None:
            pt = self.N_pt + 1
        if (x is not None) and ((pt not in self.biomarkers.index) or overwrite):
            msg = ("size of x does not match patient data structure with "
                   f"{self.N_bm} biomarkers")
            assert len(x) == self.N_bm, msg
            biomarkers.loc[pt] = x
        elif x is not None:
            warn(f"Patient {pt} exists, failing to overwrite data")
        return biomarkers

    def _append_outcomes(self, pt, tx, t, y, inplace=False, overwrite=False):
        if inplace:
            outcomes = self.outcomes
        else:
            outcomes = self.outcomes.copy()
        if ((pt, tx, t) in outcomes.index) and (not overwrite):
            warn(f"Outcome ({pt}, {tx}, {t}) exists, failing to overwrite data")
        if (pt is not None) and (tx is not None) and (t is not None):
            outcomes.loc[pt, tx, t] = y
        return outcomes
    
