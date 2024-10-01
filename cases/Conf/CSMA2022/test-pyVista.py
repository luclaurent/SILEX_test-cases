import pickle
import numpy as np

f = open('out.db','rb')

p = pickle.Unpickler(f)
data = p.load()
f.close()