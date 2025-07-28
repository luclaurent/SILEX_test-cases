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
from scipy.ndimage.filters import gaussian_filter

######################################
filename = Path(__file__).parent / 'results_parametric_tet4.pck'
f=open(filename,'rb')
frf_tet4_xfem =pickle.load(f)
f.close()

#X : lup
#Y : ldown


nbval = frf_tet4_xfem[0].shape[0]
X = frf_tet4_xfem[0].reshape((nbval, nbval))
Y = frf_tet4_xfem[1].reshape((nbval, nbval))
val_p = frf_tet4_xfem[2].reshape((nbval, nbval))
val_f = frf_tet4_xfem[3].reshape((nbval, nbval))

sigma=1.0
val_p = gaussian_filter(val_p, sigma)
sigma=1.2
val_f = gaussian_filter(val_f, sigma)


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
surf_p = ax1.plot_surface(X, Y, val_p/(1000*9.81)
                    #    , 
                    #    rstride=1, cstride=1, facecolors=rgb,
                    #    linewidth=0, antialiased=False, 
                    #    shade=False
                       )

surf_f = ax2.plot_surface(X, Y, val_f
                    #    , 
                    #    rstride=1, cstride=1, facecolors=rgb,
                    #    linewidth=0, antialiased=False, 
                    #    shade=False
                       )

# https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html
fig3, (ax3,ax4) = plt.subplots(2)
ax3.set_aspect('equal')
ax3.set_title("Elevation")
ax3.set_xlabel('up shift')
ax3.set_ylabel('down shift')
# cmap : supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 
# 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Grays', 'Grays_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 
# 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 
# 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 
# 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 
# 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 
# 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 
# 'afmhot_r', 'autumn', 'autumn_r', 'berlin', 'berlin_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 
# 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 
# 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 
# 'gist_gray_r', 'gist_grey', 'gist_grey_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 
# 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gist_yerg', 
# 'gist_yerg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'grey', 'grey_r', 'hot', 
# 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'managua', 'managua_r', 
# 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 
# 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 
# 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 
# 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'vanimo', 'vanimo_r', 
# 'viridis', 'viridis_r', 'winter', 'winter_r'

ax3.contour(X, Y, val_p, levels=10, linewidths=1)
cntr_p = ax3.contour(X, Y, val_p, levels=20, cmap="autumn_r")
fig3.colorbar(cntr_p, ax=ax3)


ax4.set_aspect('equal')
ax4.set_title("Force")
ax3.set_xlabel('up shift')
ax3.set_ylabel('down shift')
ax4.contour(X, Y, val_f, levels=5, linewidths=1)
cntr_f = ax4.contour(X, Y, val_f, levels=20, cmap="winter")
fig3.colorbar(cntr_f, ax=ax4)





plt.show()





# Z = griddata((x,y),z,(X,Y), method='cubic')

# fig = go.Figure(data=[go.Surface(z=val, x=X, y=Y)])

# fig.update_layout(fig.layout)
# fig.show()


# show_in_window(fig)



