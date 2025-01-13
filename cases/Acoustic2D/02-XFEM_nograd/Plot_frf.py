from numpy import *
import string
import time
#from scipy.sparse import *
#from scipy import *
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen
#from numpy.linalg import solve, norm
import os
#from scipy.linalg import decomp
import pylab as pl
import pickle

# Plot an other FRF for comparison 
#f=open('/home/legay/Codes/XFEM-Acoustique/V2013/test-blending/results/classic_para_results.frf','r')
#frf_ref=pickle.load(f)
#f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/classic-test1_results.frf','r')
frf_classic_2d=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test1-no-edge_results.frf','r')
frf_xfem_1_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test1-with-edge_results.frf','r')
frf_xfem_1_with_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test2-no-edge_results.frf','r')
frf_xfem_2_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test2-with-edge_results.frf','r')
frf_xfem_2_with_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test3-no-edge_results.frf','r')
frf_xfem_3_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test3-with-edge_results.frf','r')
frf_xfem_3_with_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test4-no-edge_results.frf','r')
frf_xfem_4_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test4-with-edge_results.frf','r')
frf_xfem_4_with_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test5-no-edge_results.frf','r')
frf_xfem_5_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test5-with-edge_results.frf','r')
frf_xfem_5_with_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test6-no-edge_results.frf','r')
frf_xfem_6_no_edge=pickle.load(f)
f.close()

f=open('/home/legay/Codes/XFEM-Acoustique/V2013-2D/test/results/xfem-test6-with-edge_results.frf','r')
frf_xfem_6_with_edge=pickle.load(f)
f.close()

prefsquare=20e-6*20e-6

dBerror_with_edge_1 =abs(10*log10(frf_xfem_1_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_with_edge_2 =abs(10*log10(frf_xfem_2_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_with_edge_3 =abs(10*log10(frf_xfem_3_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_with_edge_4 =abs(10*log10(frf_xfem_4_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_with_edge_5 =abs(10*log10(frf_xfem_5_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_with_edge_6 =abs(10*log10(frf_xfem_6_with_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))

globalerror_with_edge=[(sum(dBerror_with_edge_1))/400,
                       (sum(dBerror_with_edge_2))/400,
                       (sum(dBerror_with_edge_3))/400,
                       (sum(dBerror_with_edge_4))/400,
                       (sum(dBerror_with_edge_5))/400,
                       (sum(dBerror_with_edge_6))/400
                       ]
dBerror_no_edge_1 =abs(10*log10(frf_xfem_1_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_no_edge_2 =abs(10*log10(frf_xfem_2_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_no_edge_3 =abs(10*log10(frf_xfem_3_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_no_edge_4 =abs(10*log10(frf_xfem_4_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_no_edge_5 =abs(10*log10(frf_xfem_5_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))
dBerror_no_edge_6 =abs(10*log10(frf_xfem_6_no_edge[1]/prefsquare)-10*log10(frf_classic_2d[1]/prefsquare))

globalerror_no_edge=[(sum(dBerror_no_edge_1))/400,
                     (sum(dBerror_no_edge_2))/400,
                     (sum(dBerror_no_edge_3))/400,
                     (sum(dBerror_no_edge_4))/400,
                     (sum(dBerror_no_edge_5))/400,
                     (sum(dBerror_no_edge_6))/400
                     ]

nb_dll=[13,123,963,8554,77653,161196]  # nb nodes

pl.figure(1)
#pl.plot(frf_ref[0],10*log10(frf_ref[1]/prefsquare),'k-',label='Reference 3d', linewidth=3)
pl.plot(frf_classic_2d[0],10*log10(frf_classic_2d[1]/prefsquare),'k-',label='Reference', linewidth=3)

#pl.plot(frf_xfem_1_no_edge[0],10*log10(frf_xfem_1_no_edge[1]/prefsquare),'g-',label='Xfem 1 no edge', linewidth=1)
pl.plot(frf_xfem_1_with_edge[0],10*log10(frf_xfem_1_with_edge[1]/prefsquare),'b-', linewidth=1)
pl.plot(frf_xfem_1_with_edge[0][range(0,len(frf_xfem_1_with_edge[0]),20)],10*log10(frf_xfem_1_with_edge[1][range(0,len(frf_xfem_1_with_edge[1]),20)]/prefsquare),'sb', linewidth=1)
pl.plot(frf_xfem_1_with_edge[0][range(2)],10*log10(frf_xfem_1_with_edge[1][range(2)]/prefsquare),'sb-',label='Mesh 1', linewidth=1)

pl.plot(frf_xfem_2_no_edge[0],10*log10(frf_xfem_2_no_edge[1]/prefsquare),'g-', linewidth=1)
pl.plot(frf_xfem_2_no_edge[0][range(0,len(frf_xfem_2_no_edge[0]),20)],10*log10(frf_xfem_2_no_edge[1][range(0,len(frf_xfem_2_no_edge[1]),20)]/prefsquare),'<g', linewidth=1)
pl.plot(frf_xfem_2_no_edge[0][range(2)],10*log10(frf_xfem_2_no_edge[1][range(2)]/prefsquare),'<g-',label='Mesh 2, no edge enrichment', linewidth=1)

pl.plot(frf_xfem_2_with_edge[0],10*log10(frf_xfem_2_with_edge[1]/prefsquare),'m-', linewidth=1)
pl.plot(frf_xfem_2_with_edge[0][range(0,len(frf_xfem_2_with_edge[0]),20)],10*log10(frf_xfem_2_with_edge[1][range(0,len(frf_xfem_2_with_edge[1]),20)]/prefsquare),'>m', linewidth=1)
pl.plot(frf_xfem_2_with_edge[0][range(2)],10*log10(frf_xfem_2_with_edge[1][range(2)]/prefsquare),'>m-',label='Mesh 2', linewidth=1)

#pl.plot(frf_xfem_3_no_edge[0],10*log10(frf_xfem_3_no_edge[1]/prefsquare),'r-',label='Xfem 3 no edge', linewidth=1)
pl.plot(frf_xfem_3_with_edge[0],10*log10(frf_xfem_3_with_edge[1]/prefsquare),'r-', linewidth=1)
pl.plot(frf_xfem_3_with_edge[0][range(0,len(frf_xfem_3_with_edge[0]),20)],10*log10(frf_xfem_3_with_edge[1][range(0,len(frf_xfem_3_with_edge[1]),20)]/prefsquare),'or', linewidth=1)
pl.plot(frf_xfem_3_with_edge[0][range(2)],10*log10(frf_xfem_3_with_edge[1][range(2)]/prefsquare),'or-',label='Mesh 3', linewidth=1)

#pl.plot(frf_xfem_4_no_edge[0],10*log10(frf_xfem_4_no_edge[1]/prefsquare),'g-',label='Xfem 4 no edge', linewidth=1)
#pl.plot(frf_xfem_4_with_edge[0],10*log10(frf_xfem_4_with_edge[1]/prefsquare),'b-',label='Xfem 4 with edge', linewidth=1)

#pl.plot(frf_xfem_5_no_edge[0],10*log10(frf_xfem_5_no_edge[1]/prefsquare),'m-',label='Xfem 5 no edge', linewidth=1)
#pl.plot(frf_xfem_5_with_edge[0],10*log10(frf_xfem_5_with_edge[1]/prefsquare),'c-',label='Xfem 5 with edge', linewidth=1)

#pl.plot(frf_xfem_6_no_edge[0],10*log10(frf_xfem_6_no_edge[1]/prefsquare),'g-',label='Xfem 6 no edge', linewidth=1)
#pl.plot(frf_xfem_6_with_edge[0],10*log10(frf_xfem_6_with_edge[1]/prefsquare),'b-',label='Xfem 6 with edge', linewidth=1)

#pl.plot(frf_xfem_2[0],10*log10(frf_xfem_2[1]/prefsquare),'y-',label='Xfem 2', linewidth=1)
#pl.plot(frf_xfem_3[0],10*log10(frf_xfem_3[1]/prefsquare),'g-',label='Xfem 3', linewidth=1)
#pl.plot(frf_xfem_4[0],10*log10(frf_xfem_4[1]/prefsquare),'b-',label='Xfem 4', linewidth=1)
#pl.plot(frf_xfem_5[0],10*log10(frf_xfem_5[1]/prefsquare),'m-',label='Xfem 5', linewidth=1)
#pl.plot(frf_xfem_6[0],10*log10(frf_xfem_6[1]/prefsquare),'c-',label='Xfem 6', linewidth=1)

#pl.plot(frf_xfem_1_no_edge[0],10*log10(frf_xfem_1_no_edge[1]/prefsquare),'r--',label='Xfem 1 no edge', linewidth=1)
#pl.plot(frf_xfem_2_no_edge[0],10*log10(frf_xfem_2_no_edge[1]/prefsquare),'y--',label='Xfem 2 no edge', linewidth=1)
#pl.plot(frf_xfem_3_no_edge[0],10*log10(frf_xfem_3_no_edge[1]/prefsquare),'g--',label='Xfem 3 no edge', linewidth=1)
#pl.plot(frf_xfem_4_no_edge[0],10*log10(frf_xfem_4_no_edge[1]/prefsquare),'b--',label='Xfem 4 no edge', linewidth=1)
#pl.plot(frf_xfem_5_no_edge[0],10*log10(frf_xfem_5_no_edge[1]/prefsquare),'m--',label='Xfem 5 no edge', linewidth=1)
#pl.plot(frf_xfem_6_no_edge[0],10*log10(frf_xfem_6_no_edge[1]/prefsquare),'c--',label='Xfem 6 no edge', linewidth=1)

pl.axis([100.0, 300.0, 70, 130])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend(loc=3)
#'best'         : 0, (only implemented for axis legends)
#'upper right'  : 1,
#'upper left'   : 2,
#'lower left'   : 3,
#'lower right'  : 4,
#'right'        : 5,
#'center left'  : 6,
#'center right' : 7,
#'lower center' : 8,
#'upper center' : 9,
#'center'       : 10,

pl.figure(2)
pl.plot(log10(nb_dll),globalerror_with_edge,'ro--',label='With tip enrichment', linewidth=1)
pl.plot(log10(nb_dll),globalerror_no_edge,'sg--',label='With no tip enrichment', linewidth=1)
#pl.axis([1, 5, 0.0, 10])
pl.xlabel('log(nb. dofs)')
pl.ylabel('Mean dB difference')
pl.grid('on')
pl.legend()

pl.show()

