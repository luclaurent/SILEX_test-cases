import numpy as np
import string
import time
import scipy

import pylab as pl
import pickle

import sys
from pathlib import Path

import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"

import matplotlib.pyplot as plt
from matplotlib import cbook, cm
from matplotlib.colors import LightSource

######################################
filename = Path(__file__).parent / 'results_parametric_tet4.pck'
f=open(filename,'rb')
frf_tet4_xfem =pickle.load(f)
f.close()

nbval = frf_tet4_xfem[0].shape[0]
X = frf_tet4_xfem[0].reshape((nbval, nbval))
Y = frf_tet4_xfem[1].reshape((nbval, nbval))
val_p = frf_tet4_xfem[2].reshape((nbval, nbval))
val_f = frf_tet4_xfem[3].reshape((nbval, nbval))

# Set up plot
fig, (ax1,ax2) = plt.subplots(2,subplot_kw=dict(projection='3d'))

ls = LightSource(270, 45)
# To use a custom hillshading mode, override the built-in shading and pass
# in the rgb colors of the shaded surface calculated from "shade".
# rgb = ls.shade(val, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
surf = ax1.plot_surface(X, Y, val_p
                    #    , 
                    #    rstride=1, cstride=1, facecolors=rgb,
                    #    linewidth=0, antialiased=False, 
                    #    shade=False
                       )

surfb = ax2.plot_surface(X, Y, val_f
                    #    , 
                    #    rstride=1, cstride=1, facecolors=rgb,
                    #    linewidth=0, antialiased=False, 
                    #    shade=False
                       )

plt.show()





# Z = griddata((x,y),z,(X,Y), method='cubic')

# fig = go.Figure(data=[go.Surface(z=val, x=X, y=Y)])

# fig.update_layout(fig.layout)
# fig.show()


# show_in_window(fig)



