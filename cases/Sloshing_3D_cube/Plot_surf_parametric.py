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

#for i in range(len(val_f)):
#   for j in range(len(val_f[i])):
#      print(val_f[i][j])
#      if val_f[i][j]>30000:
#         val_f[i][j]=0.5*(val_f[i][j-1]+val_f[i][j+1])
#         if val_f[i][j]>30000:
#            val_f[i][j]=0.5*(val_f[i][j-1]+val_f[i][j+2])
#         print('corrigee= ',val_f[i][j])

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



# plt.show()
force_limit = 8e3
fig, ax1 = plt.subplots(1)
masked_val = np.ma.masked_where(val_f < force_limit, val_f*0.0)
ctf = ax1.contourf(X, Y, val_f )
ax1.contourf(X, Y, masked_val, levels=50,cmap="Greys", alpha=1)

# surf = ax1.imshow((val_f<force_limit).astype(int),  
#                   alpha=1,
#                   extent=(X.min(),X.max(),Y.min(),Y.max()),
#                   origin="upper", 
#                   cmap="Greys", 
#                   aspect='auto', 
#                   interpolation = 'hanning')
plt.show()

# Z = griddata((x,y),z,(X,Y), method='cubic')

# fig = go.Figure(data=[go.Surface(z=val, x=X, y=Y)])

# fig.update_layout(fig.layout)
# fig.show()


# show_in_window(fig)



