#############################################################################
#      Import libraries
#############################################################################
import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
import sys
import pylab as pl
import pickle

#############################################################################
f=open('Results_modal_pgd_frf','rb')
sol_pgd,nodes_w_pgd=pickle.load(f)
f.close()

f=open('Results_modal_classic_frf','rb')
sol_ref,nodes_w_ref=pickle.load(f)
f.close()

pl.figure(3)
pl.plot(nodes_w_pgd/(2.0*scipy.pi),scipy.log(sol_pgd),label='pgd', linewidth=2)
pl.plot(nodes_w_ref/(2.0*scipy.pi),scipy.log(sol_ref),label='ref.', linewidth=2)
pl.xlabel('freq.')
pl.ylabel('log(dep.)')
pl.grid('on')
pl.legend(loc=1)

pl.show()
