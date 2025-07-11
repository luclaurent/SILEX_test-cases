from numpy import *
import string
import time
from scipy.sparse import *
from scipy import *
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve,use_solver,minres,eigen
from numpy.linalg import solve, norm
import os
from scipy.linalg import decomp
import pylab as pl
import pickle

# Plot an other FRF for comparison 
f=open('cubeplan-classic_results.frf','r')
frf_classic_finemesh=pickle.load(f)
f.close()
f=open('cubeplan-xfem_results.frf','r')
frf_xfem_finemesh=pickle.load(f)
f.close()


#############################
#                                    Classic
#                             Coarse    Middle       Fine
# nnodes for fluid       =      2481     11674      31984
# nelem  for fluid       =     11764     62629     180971
# nnodes for structure   =       108       313        599
# nelem  for structure   =       178       560       1106
#############################
#                                      Xfem
#                              Coarse   Middle       Fine
# nnodes for fluid       =      2521     12095      33091
# nelem  for fluid       =     12636     67057     191390
# nnodes for structure   =       110       318        596
# nelem  for structure   =       182       570       1100
# nnodes fluid enriched  =       295       783       1529
# nelem  fluid enriched  =       725      2097       4223

#######################################
# Analytical frequencies of cavity 1  #
#######################################
celerity=340.0
lx=0.75
ly=0.6
lz=0.4
n=5
fluid1_freq=[]
I=[]
J=[]
K=[]
for k in range(n):
    for j in range(n):
        for i in range(n):
            fluid1_freq.append(celerity*(sqrt(((i)/lx)**2+((j)/ly)**2+((k)/lz)**2))/2)
            I.append(i)
            J.append(j)
            K.append(k)
            
#######################################
## Analytical frequencies of cavity 2 #
#######################################
#celerity=340.0
lx=0.25
ly=0.6
lz=0.4
n=5
fluid2_freq=[]
I=[]
J=[]
K=[]
for k in range(n):
    for j in range(n):
        for i in range(n):
            fluid2_freq.append(celerity*(sqrt(((i)/lx)**2+((j)/ly)**2+((k)/lz)**2))/2)
            I.append(i)
            J.append(j)
            K.append(k)

#fluid1_freq=[   226.7,   453.3,     680,   906.7,   283.3,   362.8,   534.6,   736.7,   949.9,   566.7,   610.3,   725.7,   885.2,    1069,     850,   879.7,   963.3,    1089,    1243,    1133,    1156,    1221,    1322,    1451,     425,   481.7,   621.4,   801.9,    1001,   510.8,   558.8,   682.9,   850.5,    1041,   708.3,   743.7,     841,   981.9,    1151,   950.3,     977,    1053,    1169,    1313,    1210,    1231,    1293,    1388,    1512,     850,   879.7,   963.3,    1089,    1243,     896,   924.2,    1004,    1125,    1275,    1022,    1046,    1118,    1227,    1366,    1202,    1223,    1285,    1381,    1506,    1417,    1435,    1487,    1571,    1682,    1275,    1295,    1353,    1445,    1565,    1306,    1326,    1383,    1473,    1590,    1395,    1414,    1467,    1552,    1664,    1532,    1549,    1598,    1676,    1780,    1706,    1721,    1765,    1836,    1932,    1700,    1715,    1759,    1831,    1927,    1723,    1738,    1782,    1853,    1947,    1792,    1806,    1848,    1917,    2008,    1901,    1914,    1954,    2019,    2106,    2043,    2056,    2093,    2153,    2235]
#fluid2_freq=[       680,      1360,      2040,      2720,     283.3,     736.7,      1389,      2060,      2735,     566.7,     885.2,      1473,      2117,      2778,	850,      1089,      1604,      2210,      2850,      1133,      1322,      1770,      2334,      2947,	425,     801.9,      1425,      2084,      2753,     510.8,     850.5,      1453,      2103,      2768,     708.3,     981.9,      1533,      2159,      2811,     950.3,      1169,      1659,      2250,      2881,      1210,      1388,      1821,      2372,      2977,	850,      1089,      1604,      2210,      2850,	896,      1125,      1629,      2228,      2864,      1022,      1227,      1701,      2281,      2906,      1202,      1381,      1815,      2368,      2974,      1417,      1571,      1964,      2484,      3067,      1275,      1445,      1864,      2406,      3004,      1306,      1473,      1886,      2422,      3017,      1395,      1552,      1948,      2472,      3057,      1532,      1676,      2049,      2551,      3122,      1706,      1836,      2182,      2659,      3211,      1700,      1831,      2177,      2655,      3208,      1723,      1853,      2195,      2671,      3220,      1792,      1917,      2250,      2715,      3257,      1901,      2019,      2337,      2788,      3318,      2043,      2153,      2454,      2887,      3402]

###############################
# structure eigen frequencies #
###############################
# Analytical frequencies of a structure: SS-SS-SS-SS
###############################
E=70000.0e6
nu=0.27
rho=2700.0
h=4.0e-3

D=E*h**3/(12*(1-nu**2))

n=5

strc_freq=[]
I=[]
J=[]

for j in range(n):
    for i in range(n):
        omega=sqrt( (D/(rho*h)))*( ((i+1)*pi/ly)**2+((j+1)*pi/lz)**2 )
        strc_freq.append(omega/(2*pi))
        I.append(i+1)
        J.append(j+1)

# structure numerical results

#mesh 10
#nnodes for structure= 82
strc_freq10=[   88.48476378,   174.48491316,   283.95038927,   323.7198469 ,   378.06234073,   535.42288888,   551.38965138,   643.2702044 ,   755.10542909,   784.81705762,   855.71318963,   919.79622353,  1100.09976894,  1161.85164951,  1208.976785  ,  1255.27681243,  1377.98841592,  1480.55006789,  1504.11006306,  1516.30022169,  1774.08333804,  1805.10667225,  1935.90696538,  1945.02938793,  2045.16192311,  2061.3401107 ,  2131.58360676,  2262.25319768,  2362.22598783,  2385.89707631,  2475.61497974,  2586.2887904 ,  2605.16625483,  2788.05677209,  2843.44312643,  2927.17260582,  2995.30389419,  3101.09238068,  3162.75602752,  3176.78411505]

#mesh 20
#nnodes for structure= 314
strc_freq20=[   87.02877572,   168.24473266,   270.60667264,   305.40357978,   353.39373716,   493.4646049 ,   501.01035108,   584.35195403,   669.39641119,   692.93955285,   758.10311241,   814.70380939,   956.05038137,  1018.65453345,  1040.47636476,  1085.36011163,  1129.08459458,  1289.69458418,  1294.50294946,  1281.7564811 ,  1482.63531253,  1489.91336612,  1635.35311975,  1657.04390139,  1697.71749921,  1746.13079122,  1780.65374928,  1913.93408469,  1959.24183789,  2044.87766315,  2128.89089441,  2140.42961715,  2177.57765993,  2424.69938401,  2464.1394081 ,  2518.7963905 ,  2547.33723991,  2557.98328967,  2570.18122502,  2725.08651898]

#mesh 40
#nnodes for structure= 1244
strc_freq40=[   86.69667751,   166.90972244,   267.46743218,   300.99139208,   348.06587669,   482.79053947,   489.50594053,   570.53487063,   651.8515558 ,   672.31321396,   733.56038124,   787.52607909,   917.87685926,   978.4535412 ,   999.70384275,  1033.85818124,  1081.76540619,  1226.27571419,  1218.76864909,  1219.54934476,  1392.45075829,  1412.31152969,  1530.20396099,  1557.78063453,  1580.15324987,  1641.70992347,  1664.00174974,  1779.34783117,  1809.31789085,  1894.06539258,  1976.29096979,  1970.05283472,  1998.87016602,  2232.00112825,  2251.4564764 ,  2286.36783919,  2315.62921416,  2340.35204409,  2336.47678384,  2475.38051948]

#mesh 80
#nnodes for structure= 4998
strc_freq80=[   86.61715389,   166.6119587 ,   266.68168495,   300.03314233,   346.77180835,   480.35400133,   487.02642571,   567.29993057,   647.53485989,   667.58812864,   727.76490016,   781.39830738,   908.64877583,   969.01024121,   989.22076409,  1022.45226266,  1069.67056139,  1210.49349025,  1203.72416128,  1203.84801517,  1371.5632765 ,  1392.0542371 ,  1505.98803941,  1533.47227239,  1553.48294681, 1614.10028215,  1634.2993987 ,  1748.79871761,  1775.46614201,  1856.41075799,  1937.78234303,  1930.58622251,  1957.94354449,  2180.81190309,  2201.42972264,  2234.56372746,  2261.64189223,  2282.01101875,  2282.44111422,  2417.47955261]

#mesh 160
#nnodes for structure= 19819
strc_freq160=[   86.59761554,   166.54452064,   266.49546272,   299.81375963,   346.46540015,   479.77492778,   486.43923604,   566.4483718 ,   646.46037302,   666.45986134,   726.46768271,   779.83661698,   906.5617566 ,   986.64426556,   966.60740896,  1066.6981081 ,  1019.9739571 ,  1206.82512329,  1200.17610844,  1200.15521467,  1367.04729897,  1387.08369276,  1500.56325912,  1527.34504081,  1547.32066561,  1607.47489276,  1627.47815938,  1741.05413962,  1767.73211098,  1847.89701312,  1921.43666037,  1948.14235764,  1928.13526745,  2168.74386572,  2188.84863498,  2222.1633164 ,  2248.95987127,  2268.99573813,  2269.07645657,  2402.7111029 ]

#########
# plots #
#########

fluid1_freq_plot=[[226.7,226.7,283.3,283.3,283.3,362.8,362.8,362.8,425 ,425,425 ,453.3,453.3,453.3,481.7,481.7,481.7,510.8,510.8,510.8,534.6,534.6,534.6,558.8,558.8,558.8,566.7,566.7,566.7],
                  [1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6,1e6,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6]]

fluid2_freq_plot=[[283.3,283.3,425  ,425  ,425  ,510.8,510.8,510.8,566.7,566.7],
                  [1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6 ]]

strc_freq_plot=[[86.6,86.6,166.5,166.5,166.5,266.4,266.4,266.4,299.7,299.7,299.7,346.4,346.4,346.4,479.6,479.6,479.6,486.2,486.2,486.2,566.2,566.2,566.2,646.1,646.1,646.1,666.1,666.1,666.1,726.0,726.0,726.0,779.3,779.3,779.3,905.9,905.9,905.9,965.8,965.8],
                [1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6  ,1e-6 ,1e-6 ,1e6 ]]


frf_xfem_finemesh[1][916]=0.5*(frf_xfem_finemesh[1][915]+frf_xfem_finemesh[1][917])
frf_xfem_finemesh[1][489]=0.5*(frf_xfem_finemesh[1][488]+frf_xfem_finemesh[1][490])
frf_xfem_finemesh[1][657]=0.5*(frf_xfem_finemesh[1][656]+frf_xfem_finemesh[1][658])

prefsquare=20e-6*20e-6

error=abs(10*log10(frf_classic_finemesh[1]/prefsquare)-10*log10(frf_xfem_finemesh[1]/prefsquare))

pl.figure(1)
#pl.plot(frf_classic_finemesh[0],10*log10(frf_classic_finemesh[1]/prefsquare),'k--',label='Compatible mesh', linewidth=2)
pl.plot(frf_classic_finemesh[0],10*log10(frf_classic_finemesh[1]/prefsquare),label='Compatible mesh', linewidth=2)
#pl.plot(frf_classic[0],log10(frf_classic[1]),'x--',label='Compatible mesh')
#pl.plot(frf_xfem[0],log10(frf_xfem[1]),'x-',label='Xfem approach')
pl.plot(frf_xfem_finemesh[0],10*log10(frf_xfem_finemesh[1]/prefsquare),'k-',label='XFEM approach', linewidth=2)
#pl.plot(frf_xfem_finemesh[0],10*log10(frf_xfem_finemesh[1]/prefsquare),label='XFEM approach', linewidth=2)
#pl.plot(fluid1_freq_plot[0],fluid1_freq_plot[1],'k:',label='large cavity eigen freq.')
#pl.plot(fluid1_freq_plot[0],fluid1_freq_plot[1],'--',label='large cavity eigen freq.')
#pl.plot(strc_freq120,log10(100*ones(len(strc_freq120))),'o',label='structure eigen frequencies 120')
#pl.plot(fluid2_freq_plot[0],fluid2_freq_plot[1],'k--',label='small cavity eigen freq.')
#pl.plot(fluid2_freq_plot[0],fluid2_freq_plot[1],'--',label='small cavity eigen freq.')
#pl.plot(strc_freq_plot[0],strc_freq_plot[1],'k-.',label='structure eigen freq.')
#pl.plot(strc_freq_plot[0],strc_freq_plot[1],'--',label='structure eigen freq.')
#pl.plot(fluid1_freq,log10(100*ones(len(fluid1_freq))),'o',label='large cavity eigen frequencies')
#pl.plot(fluid2_freq,log10(100*ones(len(fluid2_freq))),'o',label='small cavity eigen frequencies')
pl.axis([0.0, 600.0, 80, 150])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Mean quadratic pressure (dB)')
pl.grid('on')
pl.legend()

pl.figure(2)
pl.plot(frf_classic_finemesh[0],(error),'k-', linewidth=2)
pl.axis([0.0, 600.0, 0, 5])
pl.xlabel('Frequency (Hz)')
pl.ylabel('Error on mean quadratic pressure (dB)')
pl.grid('on')
pl.legend()

#pl.figure(2)
#pl.rc('text', usetex=True)
#pl.plot(sort(strc_freq160),'o',label='160 - 19819 nodes')
#pl.plot(sort(strc_freq80),'o',label='80 - 4998 nodes')
#pl.plot(sort(strc_freq40),'o',label='40x40 - 1244 nodes')
#pl.plot(sort(strc_freq20),'o',label='20x20 - 314 nodes')
#pl.plot(sort(strc_freq10),'o',label='10x10 - 82 nodes')
#pl.plot(sort(strc_freq),'o',label='analytic')

#pl.axis([0.0, 7.5, 0, 600])
#pl.grid('on')
#pl.legend(loc='upper left')
#pl.xlabel(r'Mode number')
#pl.savefig('tex_demo')
pl.show()

# rotation: 290 / 0 / 30
#translation: -0.1 -0.1 -0.1
# scale: 1.2 1.2 1.2
# structure: transfo, Z 0.5
# step 166 / 91.7 Hz
# step 312 / 164  Hz
# step 443 / 228  Hz
# step 515 / 263  Hz
# step 562 / 286  Hz
# step 586 / 298  Hz
