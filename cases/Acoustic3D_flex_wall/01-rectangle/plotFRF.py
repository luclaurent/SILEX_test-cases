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
for file in glob.glob("solvFRF*struct.pkl"):
    with open(file, "rb") as f:
        data = pickle.load(f)
    dataplot.append(data)
    #
    d["freq"] = np.array(data["freq"]).flatten()
    nbelem = re.findall(r'\d+', file)
    d['pquad_'+str(int(data["nelem"]))+'_'+str(int(data["nbnodes"]))] = np.array(data["mag"])
    d['Lp_'+str(int(data["nelem"]))+'_'+str(int(data["nbnodes"]))] = 10*np.log10(np.abs(np.array(data["mag"]))/pref**2)



# compute eigen frequencies
lx_large = 0.75
lx_small = 0.25
ly = 0.6
lz = 0.4
clty = 340.0
fmax = 1000
Y = 70000.0e6  # E Young
nu = 0.27  # nu
t = 4.0e-3  # thickness
rhos = 2700.0  # rho
D = (Y*t**3)/(12*(1-nu**2))   # flexural rigidity
eigfreqsmall = list()
eigfreqlarge = list()
eigfreqstruct = list()
for i in range(20):
    for j in range(20):
        for k in range(20):
            eigfreqsmall.append((clty/2)*(i**2/lx_small**2+j**2/ly**2+k**2/lz**2)**0.5)
            eigfreqlarge.append((clty/2)*(i**2/lx_large**2+j**2/ly**2+k**2/lz**2)**0.5)
            eigfreqstruct.append(np.sqrt(D/(rhos*t))*np.pi**2*((j+1)**2/ly**2+(k+1)**2/lz**2)/(2*np.pi))

eigfreqsmall = np.unique(np.array(eigfreqsmall))
eigfreqsmall = eigfreqsmall[1:50]

eigfreqlarge = np.unique(np.array(eigfreqlarge))
eigfreqlarge = eigfreqlarge[1:50]

eigfreqstruct = np.unique(np.array(eigfreqstruct))
eigfreqstruct = eigfreqstruct[0:50]


with open('3d_rectangles_eigenfreqs.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(eigfreqsmall)
    wr.writerow(eigfreqlarge)
    wr.writerow(eigfreqstruct)




df = pd.DataFrame(data=d)
df.to_csv("3d_rectangle_flex_wall_FRF.csv")

# Plot the data
fig, ax = plt.subplots()
for data in dataplot:    
    lp = 10*np.log10(np.abs(np.array(data["mag"]).flatten())/pref**2)
    ax.plot(np.array(data["freq"]).flatten(), lp, label=data["nbnodes"])
ax.legend()

for f in eigfreqsmall:
    ax.axvline(f, color='red', linestyle='--')

for f in eigfreqlarge:
    ax.axvline(f, color='blue', linestyle='--')

for f in eigfreqstruct:
    ax.axvline(f, color='green', linestyle='--')

ax.set_xlim(10, fmax)
plt.show(block=True)
plt.pause(10)
