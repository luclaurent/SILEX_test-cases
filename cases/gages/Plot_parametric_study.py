import pickle
import pylab
import scipy
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

##############################################

##ResultsFileName='Results_capteur-tri6'
##f=open(ResultsFileName+'_epsilon','rb')
##[A6,B6,E6]=pickle.load(f)
##f.close()
##
##A6,B6=scipy.meshgrid(A6,B6)
##E6=scipy.array(E6)
##
##fig=pylab.figure(1)
##ax = Axes3D(fig)
##ax.plot_wireframe(A6, B6, E6.T, rstride=1, cstride=1)
##ax.plot_surface(A6, B6 , E6.T, rstride=1, cstride=1,  alpha=0.6)
##ax.set_xlabel('a (mm)')
##ax.set_ylabel('b (mm)')
##ax.set_zlabel('Deformation mesuree par la jauge TRI6')

ResultsFileName='Results_capteur-tri3'
f=open(ResultsFileName+'_epsilon','rb')
[A3,B3,E3]=pickle.load(f)
f.close()

A3,B3=scipy.meshgrid(A3,B3)
E3=scipy.array(E3)

fig=pylab.figure(2)
ax = Axes3D(fig)
ax.plot_wireframe(A3, B3, E3.T, rstride=1, cstride=1)
ax.plot_surface(A3, B3 , E3.T, rstride=1, cstride=1,  alpha=0.6)
ax.set_xlabel('a (mm)')
ax.set_ylabel('b (mm)')
ax.set_zlabel('Deformation mesuree par la jauge TRI3')


fig=pylab.figure(3)
ax = Axes3D(fig)
##ax.plot_wireframe(A6, B6, E6.T, rstride=1, cstride=1 )
##ax.plot_surface(A6, B6 , E6.T, rstride=1, cstride=1,  alpha=1.0,color='b' ,label='tri6')
ax.plot_wireframe(A3, B3, E3.T, rstride=1, cstride=1 )
#ax.plot_surface(A3, B3 , E3.T, rstride=1, cstride=1,  alpha=0.4, cmap=cm.coolwarm)
ax.plot_surface(A3, B3 , E3.T, rstride=1, cstride=1,  alpha=1.0, color='r' ,label='tri3')
ax.set_xlabel('a (mm)')
ax.set_ylabel('b (mm)')
ax.legend()
ax.set_zlabel('Deformation mesuree par la jauge')


pylab.show()
