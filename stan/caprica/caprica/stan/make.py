import sys
from os import path
import pystan
import pickle

from caprica.models import standir

if len(sys.argv) == 2:
    targets = [sys.argv[1]]
elif len(sys.argv) > 2:
    targets = sys.argv[1:]
else:
    targets = ["lm_predict", "lm_prior", "lm_infer",
               "lf_predict", "lf_prior", "lf_infer"]
    
for prefix in targets:
    print("compiling", f"{prefix}.stan")
    try:
        sm = pystan.StanModel(path.join(standir, prefix + ".stan"))
    except ValueError as err:
        print(err)
        print(prefix, "failed to compile")
    pickle.dump(sm, open(path.join(standir, prefix + ".pkl"), "wb"))
                
