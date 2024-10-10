import scipy
import pylab
from matplotlib import pyplot

# EXPERIMENTAL : force en A
g200 = 200e-3*9.81
g100 = 100e-3*9.81
forceexpA=[g200,g200,g200]
depBexpA =[-1.58e-3,-1.50e-3,-1.64e-3]
depCexpA =[-2.34e-3,-2.26e-3,-2.28e-3]

# EXPERIMENTAL : force en B
forceexpB=[g100,g200,g200,g200]
depAexpB =[-0.9e-3,-1.65e-3,-1.50e-3,-1.64e-3]
depCexpB =[0.5e-3,0.85e-3,1.05e-3,0.95e-3]

# MODELE POUTRE, force en A
forcepoutreA = [0.0 , g200]
deppoutreA   = [0.0 , -g200*0.5**3/(3.0*70000e6*(1.0e-3*20.0e-3**3/12.0))]

# E.F. Force en A
forceefA = [ 0.0 , g200 ]
depBefA  = scipy.array([0.0, -1.51294562e-03 ])
depCefA  = scipy.array([0.0, -1.84752857e-03 ])

# E.F. Force en B
forceefB = [ 0.0 , g200 ]
depAefB  = scipy.array([0.0, -7.56472808e-04*2])
depCefB  = scipy.array([0.0, 1.43687408e-03*2])

#pointsBef = [[0.0, 0.0], [6.43609501e-01, g50], [9.68776133e-01, g50] , [0.0,0.0]]
#polygon = pyplot.Polygon()

pylab.figure(1)
pylab.plot(depBexpA,forceexpA,'ob',label='point B, experience')
pylab.errorbar(depBexpA,forceexpA,xerr=0.1e-3,yerr=1e-3*9.81, fmt='ob')
pylab.plot(depCexpA,forceexpA,'og',label='point C, experience')
pylab.errorbar(depCexpA,forceexpA,xerr=0.1e-3,yerr=1e-3*9.81, fmt='og')
pylab.plot(depBefA,forceefA,'-b',label='point B, E.F.')
pylab.plot(depCefA,forceefA,'-g',label='point C, E.F.')
pylab.plot(deppoutreA,forcepoutreA,'-r',label='modele poutre')
#pyplot.Polygon(pointsBef, fill=True)
#pylab.axis([0.0, None, 0.0, None])
pylab.xlabel('Deplacement (m)')
pylab.ylabel('Force (N)')
pylab.title('Aile en aluminium, force en A')
pylab.legend(loc=3)

pylab.figure(2)
pylab.plot(depCexpB,forceexpB,'ob',label='point C, experience')
pylab.errorbar(depCexpB,forceexpB,xerr=0.1e-3,yerr=1e-3*9.81, fmt='ob')
pylab.plot(depAexpB,forceexpB,'og',label='point A, experience')
pylab.errorbar(depAexpB,forceexpB,xerr=0.1e-3,yerr=1e-3*9.81, fmt='og')
pylab.plot(depCefB,forceefB,'-b',label='point C, E.F.')
pylab.plot(depAefB,forceefB,'-g',label='point A, E.F.')
#pylab.plot(deppoutreA,forcepoutreA,'-r',label='modele poutre')
#pyplot.Polygon(pointsBef, fill=True)
#pylab.axis([0.0, None, 0.0, None])
pylab.xlabel('Deplacement (m)')
pylab.ylabel('Force (N)')
pylab.title('Aile en aluminium, force en B')
pylab.legend(loc=3)


pylab.show()
