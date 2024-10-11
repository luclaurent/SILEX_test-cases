import scipy
import pylab
h_tri3 = scipy.array([15.0,10.0,5.0,2.0,1.0,0.5,0.3,0.2])
error_tri3 = scipy.array([0.129858538818,
                     0.0960520190949,
                     0.0510693127416,
                     0.0197007950193,
                     0.00993261608396,
                     0.00488621836765,
                     0.00293404989031,
                     0.00194947143577])

nbnodes_tri3=scipy.array([110,237,934,5973,23778,95245,264832,595775])
sig_yy_smooth_tri3=scipy.array([30.1,34,39,41.9,42.7,43.3,43.6,43.6])

#nbnodes_tet4=scipy.array([537,2567,20941,40883])
h_tet4 = scipy.array([10.0,5.0,2.0,1.0])
error_tet4 = scipy.array([0.097111090569,0.0534728316945,0.0208490158849,0.0104672015608])
sig_yy_smooth_tet4=scipy.array([38.4,40.1,43.3,43.8])

#nbnodes_tet10=scipy.array([3145,16635,])
h_tet10 = scipy.array([10.0,5.0,2.0])
error_tet10 = scipy.array([0.0313282600215,0.0121070771312,0.0039289195679])
sig_yy_smooth_tet10=scipy.array([41.5,42.4,43.3])

h_cub8 = scipy.array([10.0,5.0,2.0,1.0])
error_cub8 = scipy.array([0.0533895629703,0.0273552686629,0.010942725697,0.00545002522546])
sig_yy_smooth_cub8=scipy.array([40.4,42,43,43.4])

h_tri6= scipy.array([15,10,5,2,1,0.5])
error_tri6 = scipy.array([0.0593104600664,0.0390222956975,0.014349672345,0.0039735526227,0.00153020602994,0.000630073223845])
sig_yy_smooth_tri6 = scipy.array([33.67,36.75964292,39.81847943,42.42787186,43.26260978,43.47422746])

h_quad4 = scipy.array([15.0,10.0,5.0,2.0,1.0,0.5,0.2])*2.0
error_quad4 = scipy.array([0.130252453372,0.0931474281295,0.0511647892748,0.0209694187995,0.0104544696045,0.00520176420919,0.00207165684466])
sig_yy_smooth_quad4=scipy.array([31.53351031,36.23343943,39.98988512,42.23628186,42.98234604,43.35536588,43.57927504])

pylab.figure(1)
pylab.plot(h_tri3,error_tri3,color='b')
pylab.scatter(h_tri3,error_tri3,color='b',label='tri3')
pylab.plot(h_tet4,error_tet4,color='g')
pylab.scatter(h_tet4,error_tet4,color='g',label='tet4')
pylab.plot(h_tet10,error_tet10,color='r')
pylab.scatter(h_tet10,error_tet10,color='r',label='tet10')
pylab.plot(h_cub8,error_cub8,color='y')
pylab.scatter(h_cub8,error_cub8,color='y',label='cub8')
pylab.plot(h_tri6,error_tri6,color='c')
pylab.scatter(h_tri6,error_tri6,color='c',label='tri6')
pylab.plot(h_quad4,error_quad4,color='m')
pylab.scatter(h_quad4,error_quad4,color='m',label='quad4')
pylab.xlabel('Longueur moyenne des elements (mm)')
pylab.ylabel('Erreur globale relative')
pylab.title('Plaque trouee')
pylab.grid('on')
#pylab.axis([-1.0, 2.0, -3.0, 0.0])
pylab.yscale('log')
pylab.xscale('log')
pylab.legend()

pylab.figure(2)
pylab.plot(h_tri3,sig_yy_smooth_tri3,color='b')
pylab.scatter(h_tri3,sig_yy_smooth_tri3,color='b',label='tri3')
pylab.plot(h_tet4,sig_yy_smooth_tet4,color='g')
pylab.scatter(h_tet4,sig_yy_smooth_tet4,color='g',label='tet4')
pylab.plot(h_tet10,sig_yy_smooth_tet10,color='r')
pylab.scatter(h_tet10,sig_yy_smooth_tet10,color='r',label='tet10')
pylab.plot(h_cub8,sig_yy_smooth_cub8,color='y')
pylab.scatter(h_cub8,sig_yy_smooth_cub8,color='y',label='cub8')
pylab.plot(h_tri6,sig_yy_smooth_tri6,color='c')
pylab.scatter(h_tri6,sig_yy_smooth_tri6,color='c',label='tri6')
pylab.plot(h_quad4,sig_yy_smooth_quad4,color='m')
pylab.scatter(h_quad4,sig_yy_smooth_quad4,color='m',label='quad4')
pylab.xlabel('Longueur moyenne des elements (mm)')
pylab.ylabel('Contrainte max. sig_yy (MPa)')
pylab.title('Plaque trouee')
pylab.grid('on')
pylab.legend()

#pylab.axis([0,600000, 0.0, 44])
#pylab.yscale('log')
#pylab.xscale('log')

pylab.show()
