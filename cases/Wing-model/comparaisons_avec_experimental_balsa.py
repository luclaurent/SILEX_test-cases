import scipy
import pylab
from matplotlib import pyplot

# EXPERIMENTAL : force en A
g20 = 20e-3*9.81
g50 = 50e-3*9.81
forceexpA=[g20,g20,g20,g50,g50,g50]
depBexpA =[0.57,0.65,0.47,1.24,1.35,1.30]
depCexpA =[0.83,1.12,0.66,1.95,1.90,2.14]

# EXPERIMENTAL : force en B
forceexpB=[g20,g20,g20]
depCexpB =[1.38,1.48,1.36]

# E.F.: pour young = 4300 MPa et 3000 MPa / epaisseurs nominale 2mm et 4mm / force en A
##Number of nodes: 63684
##Number of elements: 121280
##Point A : deplacement MAXI =  [  5.68334475e-05   2.93441797e-05  -9.76052880e-04]
##Point B : deplacement MAXI =  [  5.52247310e-05   1.17179337e-05  -9.68776133e-04]
##Point C : deplacement MAXI =  [  5.49559604e-05   5.71848742e-06  -9.92153122e-04]
##Point A : deplacement MINI =  [  3.96512425e-05   2.04726835e-05  -6.80967125e-04]
##Point B : deplacement MINI =  [  3.85288821e-05   8.17530256e-06  -6.75890325e-04]
##Point C : deplacement MINI =  [  3.83413677e-05   3.98964238e-06  -6.92199852e-04]
##Number of nodes: 95192
##Number of elements: 182840
##Point A : deplacement MAXI =  [  5.69338616e-05   2.93936301e-05  -9.77696916e-04]
##Point B : deplacement MAXI =  [  5.53146710e-05   1.17376961e-05  -9.70371996e-04]
##Point C : deplacement MAXI =  [  5.50445229e-05   5.72854319e-06  -9.93903657e-04]
##Point A : deplacement MINI =  [  3.97212988e-05   2.05071838e-05  -6.82114128e-04]
##Point B : deplacement MINI =  [  3.85916309e-05   8.18909032e-06  -6.77003718e-04]
##Point C : deplacement MINI =  [  3.84031555e-05   3.99665804e-06  -6.93421156e-04]

# E.F.: pour young = 4300 MPa et 3000 MPa / epaisseurs maximale 2.1mm et 4.2mm / force en A
##Number of nodes: 63684
##Number of elements: 121280
##Point A : deplacement MAXI =  [  5.36347814e-05   2.79403041e-05  -9.28826472e-04]
##Point B : deplacement MAXI =  [  5.22396157e-05   1.10966892e-05  -9.22506951e-04]
##Point C : deplacement MAXI =  [  5.20056013e-05   5.38513730e-06  -9.42772894e-04]
##Point A : deplacement MINI =  [  3.74196149e-05   1.94932354e-05  -6.48018469e-04]
##Point B : deplacement MINI =  [  3.64462435e-05   7.74187619e-06  -6.43609501e-04]
##Point C : deplacement MINI =  [  3.62829777e-05   3.75707253e-06  -6.57748531e-04]

# E.F.: pour young = 4300 MPa et 3000 MPa / epaisseurs nominale 2mm et 4mm / force en B
##Point A : deplacement MAXI =  [ -1.21054343e-06   2.92235404e-05  -9.70371996e-04]
##Point B : deplacement MAXI =  [  5.89828221e-05   1.40686026e-05  -1.24083624e-03]
##Point C : deplacement MAXI =  [  6.90170078e-05   3.40375459e-06  -3.61274902e-04]
##Point A : deplacement MINI =  [ -8.44565185e-07   2.03885166e-05  -6.77003718e-04]
##Point B : deplacement MINI =  [  4.11508061e-05   9.81530417e-06  -8.65699704e-04]
##Point C : deplacement MINI =  [  4.81514008e-05   2.37471251e-06  -2.52052257e-04]

# MODELE POUTRE, force en A
forcepoutreA = [0.0 , g50]
deppoutreA   = [0.0 , g50*500**3/(3*4300*(4.0*20.0**3/12.0))]

# E.F. Force en A
forceefA = [ 0.0 , g50 , 0.0 , g50 ]
depBefA  = scipy.array([0.0, 9.68776133e-01 , 0.0 , 6.43609501e-01 ])*50.0/20.0
depCefA  = scipy.array([0.0, 9.92153122e-01 , 0.0 , 6.57748531e-01 ])*50.0/20.0

# E.F. Force en B
forceefA = [ 0.0 , g50 , 0.0 , g50 ]
depBefA  = scipy.array([0.0, 9.68776133e-01 , 0.0 , 6.43609501e-01 ])*50.0/20.0
depCefA  = scipy.array([0.0, 9.92153122e-01 , 0.0 , 6.57748531e-01 ])*50.0/20.0

pointsBef = [[0.0, 0.0], [6.43609501e-01, g50], [9.68776133e-01, g50] , [0.0,0.0]]
#polygon = pyplot.Polygon()

pylab.figure(1)
pylab.plot(depBexpA,forceexpA,'ob',label='point B, experience')
pylab.errorbar(depBexpA,forceexpA,xerr=0.1,yerr=1e-3*9.81, fmt='ob')
pylab.plot(depCexpA,forceexpA,'og',label='point C, experience')
pylab.errorbar(depCexpA,forceexpA,xerr=0.1,yerr=1e-3*9.81, fmt='og')
pylab.plot(depBefA,forceefA,'-b',label='point B, E.F.')
pylab.plot(depCefA,forceefA,'-g',label='point C, E.F.')
pylab.plot(deppoutreA,forcepoutreA,'-r',label='modele poutre')
#pyplot.Polygon(pointsBef, fill=True)
#pylab.axis([0.0, None, 0.0, None])
pylab.xlabel('Deplacement (mm)')
pylab.ylabel('Force (N)')
pylab.title('Force en A')
pylab.legend(loc=4)

pylab.figure(2)
pylab.plot(depCexpB,forceexpB,'ob',label='point C, experience')
pylab.errorbar(depCexpB,forceexpB,xerr=0.1,yerr=1e-3*9.81, fmt='ob')
pylab.xlabel('Deplacement (mm)')
pylab.ylabel('Force (N)')
pylab.title('Force en B')

pylab.show()
