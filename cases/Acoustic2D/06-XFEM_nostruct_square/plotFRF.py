import pickle
import matplotlib.pyplot as plt
import glob
import numpy as np
import csv
import pandas as pd
import re

pref = 20e-6

# Load the data
d = dict()
dataplot = list()
for file in glob.glob("solvFRF*.pkl"):
    with open(file, "rb") as f:
        data = pickle.load(f)
    dataplot.append(data)
    #
    d["freq"] = data["freq"]
    d['pquad_'+str(int(data["nelem"]))+'_'+str(int(data["nbnodes"]))] = data["mag"]
    d['Lp_'+str(int(data["nelem"]))+'_'+str(int(data["nbnodes"]))] = 10*np.log10(np.abs(data["mag"])/pref**2)



# compute eigen frequencies
lx = 0.8
ly = 0.6
clty = 340.0
fmax = 1000
eigfreq = list()
for i in range(10):
    for j in range(10):
        eigfreq.append((clty/2)*(i**2/lx**2+j**2/ly**2)**0.5)

eigfreq = np.sort(np.array(eigfreq))
eigfreq = eigfreq[(eigfreq<fmax) & (eigfreq>0)]

with open('2d_rectangle_eigenfreq.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(eigfreq)




df = pd.DataFrame(data=d)
df.to_csv("2d_rectangle_rigid_wall_FRF.csv")

# Plot the data
fig, ax = plt.subplots()
for data in dataplot:    
    lp = 10*np.log10(np.abs(data["mag"])/pref**2)
    ax.plot(data["freq"], lp, label=data["nbnodes"])
ax.legend()

for f in eigfreq:
    ax.axvline(f, color='red', linestyle='--')


plt.show(block=True)
plt.pause(10)
