import scipy
import pylab

h_tet4        =[5.0/10.0,3.0/10.0,1.0/10.0,0.6/10.0]
nb_nodes_tet4 =[1011,3065,35021,161639]
nb_elem_tet4  =[3622,11891,154526,865338]
error_tet4    =[0.304807012325,0.23422662453,0.124036045534,0.0853394148]
cpu_time_tet4 =[0.03,0.18,35.92,551.6]

coeff_cpu=551.6/495.169086

h_tet10       =[6.0/10.0,3.0/10.0,1.5/10.0]
nb_nodes_tet10=[4089,21379,132484]
nb_elem_tet10 =[2004,12339,85592]
error_tet10   =[0.22260924136651447,0.08968592685169426,0.04198840126250299]
cpu_time_tet10=[2.2989799999999994*coeff_cpu,20.867652*coeff_cpu,369.966401*coeff_cpu]

h_pente_1     =[1,0.01]
error_pente_1 =[1,0.01]
h_pente_2     =[1,0.1]
error_pente_2 =[1,0.01]

pente_erreur_tet4=(scipy.log(error_tet4[0])-scipy.log(error_tet4[2]))/(scipy.log(h_tet4[0])-scipy.log(h_tet4[2]))
pente_erreur_tet10=(scipy.log(error_tet10[0])-scipy.log(error_tet10[2]))/(scipy.log(h_tet10[0])-scipy.log(h_tet10[2]))

print('pente erreur tet4  = ',pente_erreur_tet4)
print('pente erreur tet10 = ',pente_erreur_tet10)

sig_max_tet4        =[86.8,81.5,107,126]
sig_max_tet10       =[69.3+0.1,72.9+0.1,75.9]
sig_max_lissee_tet4 =[53.9,56.9,76.7,83.3]
sig_max_lissee_tet10=[69.3,72.9,75.7]
dep_tet4            =[0.0207,0.0224,0.0241,0.0244]
dep_tet10           =[0.0244,0.0246,0.0247]

pylab.figure(1)
pylab.plot(h_tet4,error_tet4,'o-',color='b',label='tet-4')
pylab.plot(h_tet10,error_tet10,'^-',color='g',label='tet-10')
pylab.plot(h_pente_1,error_pente_1,'-',color='r',label='erreur optimale de pente 1')
pylab.plot(h_pente_2,error_pente_2,'-',color='k',label='erreur optimale de pente 2')
#pylab.axis([-1.0, 2.0, -3.0, 0.0])
pylab.yscale('log')
pylab.xscale('log')
pylab.xlabel('Taille relative des elements (h)')
pylab.ylabel('Erreur globale ZZ1')
pylab.legend(loc=4)
pylab.grid('on')



pylab.figure(2)
pylab.plot(nb_nodes_tet4,dep_tet4,'o-',color='b',label='tet-4')
pylab.plot(nb_nodes_tet10,dep_tet10,'^-',color='g',label='tet-10')
pylab.axis([1000.0, 200000.0, 0.02, 0.025])
pylab.xscale('log')
pylab.xlabel('Nb. noeuds')
pylab.ylabel('Deplacement maximal (mm)')
pylab.legend(loc=4)
pylab.grid('on')
#pylab.axis([1000, 200000, 0.0, 0.03])


pylab.figure(3)
pylab.plot(nb_nodes_tet4,sig_max_tet4,'o--',color='b',label='Contrainte dans les elements, tet-4')
pylab.plot(nb_nodes_tet4,sig_max_lissee_tet4,'o-',color='b',label='Contrainte lissee, tet-4')
pylab.plot(nb_nodes_tet10,sig_max_tet10,'^--',color='g',label='Contrainte dans les elements, tet-10')
pylab.plot(nb_nodes_tet10,sig_max_lissee_tet10,'^-',color='g',label='Contrainte lissee, tet-10')
pylab.xscale('log')
pylab.xlabel('Nb. noeuds')
pylab.ylabel('Contrainte VM (MPa)')
pylab.legend(loc=4)
pylab.grid('on')
pylab.axis([1000, 200000, 0.0, 130])

pylab.figure(4)
pylab.plot(nb_nodes_tet4,cpu_time_tet4,'o-',color='b',label='tet-4')
pylab.plot(nb_nodes_tet10,cpu_time_tet10,'^-',color='g',label='tet-10')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('Nb. noeuds')
pylab.ylabel('Temps de resolution (s)')
pylab.legend(loc=4)
pylab.grid('on')


pylab.show()
