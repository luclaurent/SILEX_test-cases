import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path

filename = Path(__file__).parent / 'results_parametric_tet4.pck'

with open(filename, "rb") as f:
    X, Y, val_p, val_f, all_frf, all_forces = pickle.load(f)

pl.figure(1)
for frf in all_frf:
    pl.plot(frf)
    
    
pl.figure(2)
for frf in all_forces:
    pl.plot(frf)
    
pl.show()
