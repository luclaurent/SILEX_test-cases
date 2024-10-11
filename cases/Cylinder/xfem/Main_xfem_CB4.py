import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg
#from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve,use_solver,minres,eigen,cg
#from numpy.linalg import solve, norm

#import os
#from scipy.linalg import decomp
import pylab as pl
import pickle
import mumps 

import sys
sys.path.append('../../../librairies')

import silex_lib_xfem_acou_tet4
import silex_lib_dkt
import silex_lib_gmsh

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0

mycomm=comm_mumps_one_proc()

# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 4  python3 Main_toto.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py
#

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################
time_init=time.ctime()
tic00=time.clock()
if rank==0:
    print ("time at the beginning of the computation:",time.ctime())

# parallepipedic cavity with plane structure
mesh_file='geom/cyl'

# number of modes

##nb_mode_F=38
##nb_mode_S=24
##results_file='results/xfem_struc_CB_projection_fluid_damping_38F-24S'

##nb_mode_F=102
##nb_mode_S=32
##results_file='results/xfem_struc_CB_projection_fluid_damping_102F-32S'

# values up to 1000Hz
#nb_mode_F=35
#nb_mode_S=17

# values up to 1000Hz
# nb_mode_S_enc_k=17
# nb_mode_S_fre_k=21

# values up to 500Hz
nb_mode_F=35
nb_mode_S=17
nb_mode_S_enc_k=6
nb_mode_S_fre_k=6

freq_struc_fixed= [3.2906175167313334, 8.180179186217071, 15.873102174779305, 50.14043244761897, 52.993238050815087, 102.56463269370788, 172.72714269115841, 190.79511610462433, 260.66391519364248, 366.34086322227813, 432.99576809098863, 489.83009059804635, 631.11794403836245, 743.0579531579715, 790.20320429023457, 967.13846133853986, 1038.9961760800286, 1160.5149852433867, 1264.7767976839625, 1357.4880400183185, 1396.9551391646446, 1434.9110000748217, 1619.4161718039863, 1655.51797371045, 1867.7066887498565, 1868.5211453942854, 2114.0139520780285, 2137.3757389778843, 2424.7885875709621, 2428.9056540966353, 2611.9023136348142, 2733.8434716351353, 3015.725955917625, 3056.9789688076539, 3184.3458489892382, 3200.9541323953804, 3411.9392545229521, 3567.4816731084129, 3776.8702301329495, 3861.7819953135208, 4115.4364344148225, 4161.5840925369912, 4568.0686295150199, 4583.5602332911412, 4689.2283802560951, 4995.2571156985177, 5166.0253950526012, 5227.5682979032135, 5387.826275644964, 5446.9258810030724, 5849.3036034851921, 5918.0381131562826, 6131.7353521919167, 6401.8033950419567, 6498.0848406737696, 6859.6816170487582, 6913.3226980752079, 7164.6068427839164, 7219.0303539756042, 7451.4013407729208, 7564.2732841931338, 8003.1436321381952, 8028.1468256363505, 8236.5338447673676, 8583.7819338909994, 8861.228771987935, 8942.1755864096231, 9163.614270933349, 9193.8288148921019, 9518.3518521787137, 9803.4301295209316, 9942.1411043198041, 10169.833579372385, 10442.310219274123, 10806.713280037111, 11017.062114907834, 11115.691291227773, 11196.27106863426, 11453.980913622461, 11804.812485668652, 12125.669517215973, 12146.984147522486, 12524.274289396308, 12831.399712221784, 13211.392373242088, 13267.769124023289, 13270.83849437748, 13507.085902313844, 14020.080540990342, 14250.408695611268, 14455.824837351875, 14808.939950057626, 14969.297841702612, 15237.289245463549, 15620.999685516215, 15660.5501107326, 15717.291569976807, 16443.457680197222, 16505.443670305303, 16904.046288652698, 0.0]
freq_fluid      = [1.8319961067894154e-05, 170.02422340026291, 170.02449860109695, 240.48522473103102, 340.1939720159765, 340.19548240935183, 380.40264706657763, 380.40418084245499, 481.38683471606276, 510.65577408937793, 510.66324079262751, 538.3529970577199, 538.36576280837892, 614.08363415088047, 614.09721377113829, 681.54150589090023, 681.56568446506878, 702.60468424671876, 702.65923221076514, 723.10161617646224, 762.43486926276125, 762.45654395562997, 853.0184302677369, 853.03862792721122, 853.05050712675961, 853.05502939888163, 870.03226997677245, 870.08425453824532, 919.29360828466088, 919.29769188394198, 966.08411549437574, 996.06742277193575, 996.09631494629252, 1025.2784673714771, 1025.3198911289346, 1039.5104733373753, 1039.640362608756, 1081.3693817240585, 1081.4273166526548, 1094.9332481441531, 1094.9624424406118, 1147.6907629212924, 1147.7642316590798, 1198.4074954141506, 1198.4330426890399, 1210.6049717026674, 1210.7793437466685, 1210.8486761739503, 1234.9453479282172, 1235.0476805158596, 1246.9695033869048, 1247.1421313458013, 1305.3507567373701, 1305.5116648121732, 1339.2186585738907, 1339.4126475768942, 1372.4924962029881, 1372.6528574061279, 1383.2730263446454, 1383.3551669622202, 1383.4424322926848, 1383.5692222729826, 1415.6837479365593, 1415.721696615258, 1457.21429796243, 1467.6297290258569, 1467.8344183839979, 1477.8754397720327, 1477.9397586058371, 1537.9079020218323, 1538.2547079715594, 1547.7360963650324, 1547.9974677994244, 1557.3494275025587, 1557.7551926551332, 1586.2855383147423, 1586.3600611856561, 1586.6633496276224, 1586.8802069242686, 1624.3100142614937, 1624.4576527688905, 1633.5995926084186, 1633.7991254278677, 1697.7493934637787, 1697.7806666840916, 1706.5743822187069, 1724.2878697133506, 1724.378303225217, 1724.613971629344, 1724.7133181342156, 1733.0689240581075, 1733.8479254707465, 1759.7688100054502, 1759.9404350689551, 1777.0444405396233, 1777.0969566010631, 1802.6504724575602, 1802.8104566579025, 1836.2446456258151, 1836.5904214581783, 1861.4280044749942, 1861.6873864912295, 1869.7708017361188, 1869.9200493497788, 1902.6888085694638, 1902.7634549590884, 1910.4877238733825, 1911.2166210448854, 1934.4109565073577, 1934.9043208965754, 1935.2409173748495, 1935.4734180259231, 1958.8554565188895, 1974.3385799518182, 1974.4736391621771, 1974.5718769886048, 1974.9598870500697, 2021.2499108721488, 2021.3935846184534, 2028.8927710852445]
freq_struc_free = [0.00042726944161095546, 0.0003337369160671701, 0.00026498904656790798, 0.00031697609363126556, 0.28117358390922309, 0.64543585990435759, 18.405277019589029, 52.490533355551783, 104.89942098663195, 143.54645128107603, 175.05139276744336, 262.942433194827, 364.63801210281565, 368.60688774177243, 492.11978975640722, 633.52166193833034, 658.32344890837987, 792.89748879650836, 912.40156863915877, 970.54051704043866, 1146.3225971091217, 1166.11855900413, 1281.7160003108943, 1380.5233331769741, 1536.3574022919674, 1612.806455061025, 1723.1160281510472, 1864.6793635562069, 2006.101900208648, 2135.3263395734207, 2228.9342832543252, 2265.1604927809999, 2423.5310519215654, 2555.0034042737047, 2733.0599317661954, 2746.6482057426606, 3059.2728888868105, 3191.1053002032786, 3269.9661769806571, 3407.2691119931969, 3773.8909598236928, 3807.4980205472752, 3888.6946272218679, 4159.1101543286923, 4176.8273624447093, 4359.3146625089057, 4566.1458823760604, 4610.7439405567984, 4956.1008536849085, 4994.8578872299167, 5296.8061139812526, 5443.0018920277762, 5638.3071991017914, 5914.3419662694159, 5936.2665837292179, 6162.1334083702441, 6386.5266088978833, 6398.9534401095671, 6574.6114421946559, 6911.1212116202414, 7126.023229242337, 7280.8487217056645, 7447.2254873462234, 7807.5150517403072, 7999.9373401098537, 8084.9481228893119, 8167.6941651425004, 8484.633276226441, 8580.3800085671774, 8963.0008314581282, 9151.3559391986819, 9178.3098473877435, 9724.6067992251574, 9799.3489980302984, 10051.121950837918, 10181.548656737958, 10388.818298963095, 10439.241923513682, 11053.452940358917, 11098.024770025524, 11114.332104784484, 11711.239204481737, 11800.068827398411, 12199.548467363347, 12212.314242375651, 12388.303708354839, 12519.655427432721, 13056.134950860722, 13260.322945471815, 13378.403911092764, 13752.345537240053, 14014.892526903986, 14220.840846889942, 14499.008100240122, 14561.39643652546, 14804.276631518838, 15242.845526668087, 15615.499822415561, 15781.563236703883, 15964.724867454284]

if rank==0:
    print("last freq / fluid / ",freq_fluid[nb_mode_F-1])
    print("last freq / struc / ",freq_struc_fixed[nb_mode_S-1])

    if nb_mode_S_enc_k!=0:
        print("last freq / interface, fixed modes / ",freq_struc_fixed[nb_mode_S_enc_k-1])
    else:
        print("last freq / interface, fixed modes / None")
        
    if nb_mode_S_fre_k!=0:
        print("last freq / interface, free  modes / ",freq_struc_free[nb_mode_S_fre_k-1])
    else:
        print("last freq / interface, free  modes / None")

nb_mode_A=nb_mode_S_enc_k+nb_mode_S_fre_k+1 

freq_shift=113

#results_file='results/xfem_struc_CB4_'+str(nb_mode_F)+'F-'+str(nb_mode_S)+'S-'+str(nb_mode_S_enc_k)+'Aenc-'+str(nb_mode_S_fre_k)+'Afree'
results_file='results/xfem_struc_CB4'

# Fluid: air
celerity=340.0
rho=1.2

# Fluid: water
##celerity=1500.0
##rho=1000.0

# shell structure
material=[]
material.append(70000.0e6)
material.append(0.27)
material.append(6.0e-3)
material.append(2700.0)

freq_ini     = 5.0
freq_end     = 1000.0
nb_freq_step_per_proc=80

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=0

# structure damping
modal_damping_S=0.02

# fluid damping
#fluid_damping=1.0
fluid_damping=(1.0+0.01j)

##############################################################
# Load fluid mesh
##############################################################
tic = time.clock()

tic0=time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_fluid.msh',3)
fluid_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',4,1)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,4)
if rank==0:
    print ("nnodes for fluid=",fluid_nnodes)
    print ("nelem for fluid=",fluid_nelem)

##############################################################
# Load structure mesh
##############################################################
struc_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_struc.msh',3)
struc_elements,tmp = silex_lib_gmsh.ReadGmshElements(mesh_file+'_struc.msh',2,2)

struc_nnodes   = struc_nodes.shape[0]
struc_nelem    = struc_elements.shape[0]
struc_ndof     = struc_nnodes*6

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_struc_mesh',struc_nodes,struc_elements,2)

if rank==0:
    print ("nnodes for structure=",struc_nnodes)
    print ("nelem for structure=",struc_nelem)


##############################################################
# Material, Boundary conditions
##############################################################
SolvedDofF=list(range(fluid_ndof))

# Find the fixed dofs and the free dofs of the structure

tmp=scipy.sparse.find(struc_nodes[:,0]==0.0)# x=0
FixedStrucNodes=tmp[1]+1
FixedStrucDofUx=(FixedStrucNodes-1)*6
FixedStrucDofUy=(FixedStrucNodes-1)*6+1
FixedStrucDofUz=(FixedStrucNodes-1)*6+2
FixedStrucDofRx=(FixedStrucNodes-1)*6+3
FixedStrucDofRy=(FixedStrucNodes-1)*6+4
FixedStrucDofRz=(FixedStrucNodes-1)*6+5

FixedStrucDof=scipy.hstack([FixedStrucDofUx,FixedStrucDofUy,FixedStrucDofUz,FixedStrucDofRx,FixedStrucDofRy,FixedStrucDofRz])

SolvedDofS=scipy.setdiff1d(range(struc_ndof),FixedStrucDof)

# To impose the load on the structure
tmp=scipy.sparse.find(struc_nodes[:,1]==0.0)# y=0.0
IdNodeLoadStructure=tmp[1]+1

FS=scipy.zeros(struc_ndof)

for i in range(len(IdNodeLoadStructure)-1):
    IddofLoadStructure=[(IdNodeLoadStructure[i]-1)*6]
    FS[IddofLoadStructure]=FS[IddofLoadStructure]+1.0/len(IdNodeLoadStructure)
    IddofLoadStructure=[(IdNodeLoadStructure[i+1]-1)*6]
    FS[IddofLoadStructure]=FS[IddofLoadStructure]+1.0/len(IdNodeLoadStructure)

##################################################################
# compute level set
##################################################################
tic = time.clock()

LevelSet,distance = silex_lib_xfem_acou_tet4.computelevelset(fluid_nodes,struc_nodes,struc_elements)


toc = time.clock()
if rank==0:
    print ("time to compute level set:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_signed_distance',fluid_nodes,fluid_elements,4,[[[LevelSet],'nodal',1,'Level set']])

##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()

LSEnrichedElements,NbLSEnrichedElements=silex_lib_xfem_acou_tet4.getenrichedelementsfromlevelset(fluid_elements,LevelSet)
LSEnrichedElements=LSEnrichedElements[list(range(NbLSEnrichedElements))]


EnrichedElements,NbEnrichedElements=silex_lib_xfem_acou_tet4.getsurfenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements[LSEnrichedElements])
EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])
EnrichedElements=LSEnrichedElements[EnrichedElements-1]
toc = time.clock()
if rank==0:
    print ("time to find surface enriched elements:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_LSenriched_elements',fluid_nodes,fluid_elements[LSEnrichedElements],4)
    silex_lib_gmsh.WriteResults2(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],4)

##############################################################
# Compute Standard Fluid Matrices
##############################################################
tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_lib_xfem_acou_tet4.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time to compute fluid matrices:",toc-tic)

##############################################################
# Compute structure matrices
##############################################################
tic = time.clock()

IIks,JJks,Vks,Vms=silex_lib_dkt.stiffnessmatrix(struc_nodes,struc_elements,material)

KSS = scipy.sparse.csc_matrix( (Vks,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )
MSS = scipy.sparse.csc_matrix( (Vms,(IIks,JJks)), shape=(struc_ndof,struc_ndof) )

toc = time.clock()
if rank==0:
    print ("time for computing structure:",toc-tic)

##################################################################
# Compute Heaviside enrichment
##################################################################
tic = time.clock()

Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_xfem_acou_tet4.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,celerity,rho)

KAAheaviside = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
MAAheaviside = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
KAFheaviside = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
MAFheaviside = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofA=Enrichednodes-1

KAA=KAAheaviside
MAA=MAAheaviside
KAF=KAFheaviside
MAF=MAFheaviside

toc = time.clock()
if rank==0:
    print ("time to compute Heaviside enrichment:",toc-tic)

##################################################################
# Compute coupling terms on interface
##################################################################
tic = time.clock()

IIc1,JJc1,Vc1=silex_lib_xfem_acou_tet4.computexfemcoupling1(fluid_nodes,struc_nodes,fluid_elements,struc_elements,EnrichedElements)
IIc2,JJc2,Vc2=silex_lib_xfem_acou_tet4.computexfemcoupling2(fluid_nodes,struc_nodes,fluid_elements,struc_elements,EnrichedElements,LevelSet)

CSA=0.5*scipy.sparse.csc_matrix( (Vc1,(IIc1,JJc1)), shape=(struc_ndof,fluid_ndof) )+0.5*scipy.sparse.csc_matrix( (Vc2,(IIc2,JJc2)), shape=(struc_ndof,fluid_ndof) )

toc = time.clock()
if rank==0:
    print ("time to compute coupling matrices:",toc-tic)

##################################################################
# Compute eigen modes of the structure
##################################################################
tic = time.clock()

nb_mode_S_computed=max(nb_mode_S,nb_mode_S_enc_k)

eigen_values_S_computed,eigen_vectors_S_computed= scipy.sparse.linalg.eigsh(KSS[SolvedDofS,:][:,SolvedDofS],nb_mode_S_computed,MSS[SolvedDofS,:][:,SolvedDofS],sigma=0,which='LM')

eigen_values_S  = eigen_values_S_computed[list(range(nb_mode_S))]
eigen_vectors_S = eigen_vectors_S_computed[:,list(range(nb_mode_S))]

toc = time.clock()
if rank==0:
    print ("time for computing the structure modal basis (fixed structure):",toc-tic)

tic = time.clock()
if nb_mode_S_fre_k!=0:
    eigen_values_free_S,eigen_vectors_free_S= scipy.sparse.linalg.eigsh(KSS,nb_mode_S_fre_k,MSS,sigma=0,which='LM')
    freq_eigv_free_S=list(scipy.sqrt(abs(eigen_values_free_S))/(2*scipy.pi))

toc = time.clock()
if rank==0:
    print ("time for computing the structure modal basis (free structure):",toc-tic)

##################################################################
# Add static solution to the structure basis
##################################################################
tic = time.clock()
Static_mode_S = mumps.spsolve( KSS[SolvedDofS,:][:,SolvedDofS] , FS[SolvedDofS] , comm=mycomm ).T
#Static_mode_S = scipy.sparse.linalg.spsolve( KSS[SolvedDofS,:][:,SolvedDofS] , FS[SolvedDofS] ).T
toc = time.clock()
if rank==0:
    print ("time for computing the structure static mode:",toc-tic)

##################################################################
# Orthogonalisation of the static mode
##################################################################
tic = time.clock()
eigen_vectors_S=scipy.sparse.csc_matrix(eigen_vectors_S)

lines=list(range(len(SolvedDofS)))
S=scipy.sparse.coo_matrix(Static_mode_S)

MSS_static_S1 = scipy.dot(MSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT__MSS_static_11 = scipy.array(S*MSS_static_S1.todense())[0][0]
S=S/(scipy.sqrt(staticT__MSS_static_11))

for i in range(nb_mode_S):
    a=eigen_vectors_S[lines,:][:,i]
    b=KSS[SolvedDofS,:][:,SolvedDofS]*S.T
    tmp = scipy.array(a.T*b.todense())[0][0]/eigen_values_S[i]
    S=S-tmp*a.T
    MSS_static_S1 = MSS[SolvedDofS,:][:,SolvedDofS]*S.T
    staticT__MSS_static_11 = scipy.array(S*MSS_static_S1.todense())[0][0]
    S=S/(scipy.sqrt(staticT__MSS_static_11))

toc = time.clock()
if rank==0:
    print ("time for orthogonalization of the static mode:",toc-tic)

##################################################################
# Build and save the Structure Basis
##################################################################
tic = time.clock()

PSn = scipy.sparse.construct.bmat( [ [scipy.sparse.csc_matrix(eigen_vectors_S),scipy.sparse.csc_matrix(S).T] ] )

freq_eigv_S=list(scipy.sqrt(eigen_values_S)/(2*scipy.pi))
freq_eigv_S.append(0.0)

eigen_vector_S_list=[]
for i in range(PSn.shape[1]):
    Q=scipy.zeros(struc_ndof)
    Q[SolvedDofS]=PSn.todense()[:,i]
    disp=scipy.zeros((struc_nnodes,3))
    disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
    disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
    disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
    eigen_vector_S_list.append(disp)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults2(results_file+'_structure_modes',struc_nodes,struc_elements,2,[[eigen_vector_S_list,'nodal',3,'modes']])


##################################################################
# Compute structure damping matrix
##################################################################
VDnn = 2.0*modal_damping_S*scipy.sqrt(eigen_values_S)
IIDnn = list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+nb_mode_S))
JJDnn = list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+nb_mode_S))

D = scipy.sparse.coo_matrix( (VDnn,(IIDnn,JJDnn)), shape=(nb_mode_F+nb_mode_A+nb_mode_S+1,nb_mode_F+nb_mode_A+nb_mode_S+1) )


##################################################################
# Compute eigen modes of the fluid
##################################################################
tic = time.clock()

eigen_values_F,eigen_vectors_F= scipy.sparse.linalg.eigsh(KFF[SolvedDofF,:][:,SolvedDofF],nb_mode_F,MFF[SolvedDofF,:][:,SolvedDofF],sigma=0,which='LM')

freq_eigv_F=list(scipy.sqrt(eigen_values_F)/(2*scipy.pi))

if (flag_write_gmsh_results==1) and (rank==0):
    eigen_vector_F_list=[]
    for i in range(nb_mode_F):
        tmp=eigen_vectors_F[:,i].real
        eigen_vector_F_list.append(tmp[SolvedDofF])
    silex_lib_gmsh.WriteResults2(results_file+'_fluid_modes',fluid_nodes,fluid_elements,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

toc = time.clock()
if rank==0:
    print ("time for computing the fluid modes:",toc-tic)


##################################################################
# Compute Psi_Fk and Psi_Ak for the fluid
##################################################################
tic = time.clock()

Psi_Fk=scipy.zeros((len(SolvedDofF),nb_mode_A))
Psi_Ak=scipy.zeros((len(SolvedDofA),nb_mode_A))
ModeS=scipy.zeros((struc_ndof,nb_mode_A))

K=scipy.sparse.construct.bmat( [[KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
                                [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]
                                ] )

M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                ] )
#nb_mode_S_enc_k
#nb_mode_S_fre_k

for i in range(nb_mode_S_enc_k):
    #ModeS[SolvedDofS,i]=scipy.array(PSn.todense()[:,i].T)[0]
    ModeS[SolvedDofS,i]=eigen_vectors_S_computed[:,i]

ModeS[SolvedDofS,nb_mode_S_enc_k]=scipy.array(S.todense())[0]
#scipy.sparse.csc_matrix(S).T

for i in range(nb_mode_S_fre_k):
    ModeS[:,i+nb_mode_S_enc_k+1]=eigen_vectors_free_S[:,i]

toc = time.clock()
if rank==0:
    print ("time to form the system to solve PSI_Fk and Psi_Ak:",toc-tic)

omega_cst=freq_shift*2.0*scipy.pi

tic = time.clock()
MySolve = scipy.sparse.linalg.factorized( scipy.sparse.csc_matrix(K-omega_cst*omega_cst*M) ) # Makes LU decomposition.
toc = time.clock()
if rank==0:
    print ("time to factorized the system to solve PSI_Fk and Psi_Ak:",toc-tic)

tic = time.clock()
for k in range(nb_mode_A):
    tmp=scipy.hstack( [scipy.zeros(len(SolvedDofF)),CSA[:,SolvedDofA].T*ModeS[:,k]] )
    Xi=MySolve( tmp)
    #Xi = mumps.spsolve( scipy.sparse.csc_matrix(K-OmegS[k]*OmegS[k]*M) , tmp , comm=mycomm )
    #Xi = mumps.spsolve( scipy.sparse.csc_matrix(K-omega_cst*omega_cst*M) , tmp , comm=mycomm )
    #if rank==0:
    #    print("k=",k)

    Psi_Fk[:,k]=Xi[list(range(len(SolvedDofF)))]
    Psi_Ak[:,k]=Xi[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]

#######
#test:

import gram_schmidt
Psi=scipy.vstack([Psi_Fk,Psi_Ak])
import numpy

Psi_ortho=gram_schmidt.GS(Psi)
Psi_Fk=Psi_ortho[list(range(len(SolvedDofF))),:]
Psi_Ak=Psi_ortho[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA))),:]

#######

#Psi_Fk=scipy.sparse.csc_matrix(Psi_Fk)
#Psi_Ak=scipy.sparse.csc_matrix(Psi_Ak)

toc = time.clock()
if rank==0:
    print ("time to compute PSI_Fk and Psi_Ak:",toc-tic)

if (flag_write_gmsh_results==1) and (rank==0):
    eigen_vector_F_list=[]
    for i in range(nb_mode_A):
        tmp=Psi_Fk[:,i].real
        tmp[SolvedDofA]=tmp[SolvedDofA]+scipy.sign(LevelSet[SolvedDofA])*Psi_Ak[:,i].real
        eigen_vector_F_list.append(tmp)
        #CorrectedPressure[SolvedDofA]=(CorrectedPressure[SolvedDofA].T+scipy.array(enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)).T

    silex_lib_gmsh.WriteResults2(results_file+'_PsiFk',fluid_nodes,fluid_elements,4,[[eigen_vector_F_list,'nodal',1,'pressure']])

##################################################################
# Construct the whole system
##################################################################

# Fluid part
tic = time.clock()

VK_diag_pp = abs(eigen_values_F)
VM_diag_pp = eigen_values_F/eigen_values_F
IIDpp = list(range(nb_mode_F))
JJDpp = list(range(nb_mode_F))

K_diag_pp= scipy.sparse.csc_matrix( (VK_diag_pp,(IIDpp,JJDpp)), shape=(nb_mode_F,nb_mode_F) )
M_diag_pp= scipy.sparse.csc_matrix( (VM_diag_pp,(IIDpp,JJDpp)), shape=(nb_mode_F,nb_mode_F) )


K_pk_1 = scipy.dot(scipy.array(eigen_vectors_F.T*KFF[SolvedDofF,:][:,SolvedDofF]),Psi_Fk)
K_pk_2 = scipy.dot(scipy.array(eigen_vectors_F.T*KAF[SolvedDofF,:][:,SolvedDofA]),Psi_Ak)
K_pk = K_pk_1+K_pk_2

K_kk_1 = scipy.dot(Psi_Fk.T*KFF[SolvedDofF,:][:,SolvedDofF],Psi_Fk)
K_kk_2 = scipy.dot(Psi_Fk.T*KAF[SolvedDofF,:][:,SolvedDofA],Psi_Ak)
K_kk_3 = K_kk_2.T
K_kk_4 = scipy.dot(Psi_Ak.T*KAA[SolvedDofA,:][:,SolvedDofA],Psi_Ak)

K_kk = K_kk_1+K_kk_2+K_kk_3+K_kk_4

CnA = PSn.T*CSA[SolvedDofS,:][:,SolvedDofA]
Cnk = CnA*Psi_Ak

toc = time.clock()
if rank==0:
    print ("time to compute K projections:",toc-tic)

tic = time.clock()

M_pk_1 = scipy.dot(scipy.array(eigen_vectors_F.T*MFF[SolvedDofF,:][:,SolvedDofF]),Psi_Fk)
M_pk_2 = scipy.dot(scipy.array(eigen_vectors_F.T*MAF[SolvedDofF,:][:,SolvedDofA]),Psi_Ak)
M_pk = M_pk_1+M_pk_2

M_kk_1 = scipy.dot(Psi_Fk.T*MFF[SolvedDofF,:][:,SolvedDofF],Psi_Fk)
M_kk_2 = scipy.dot(Psi_Fk.T*MAF[SolvedDofF,:][:,SolvedDofA],Psi_Ak)
M_kk_3 = M_kk_2.T
M_kk_4 = scipy.dot(Psi_Ak.T*MAA[SolvedDofA,:][:,SolvedDofA],Psi_Ak)

M_kk = M_kk_1+M_kk_2+M_kk_3+M_kk_4

toc = time.clock()
if rank==0:
    print ("time to compute M projections:",toc-tic)

tic = time.clock()

eigen_vectors_F=scipy.sparse.csc_matrix(eigen_vectors_F)

VK_diag_nn = eigen_values_S
VM_diag_nn = eigen_values_S/eigen_values_S
IIDnn = list(range(nb_mode_S))
JJDnn = list(range(nb_mode_S))

K_diag_nn= scipy.sparse.csc_matrix( (VK_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )
M_diag_nn= scipy.sparse.csc_matrix( (VM_diag_nn,(IIDnn,JJDnn)), shape=(nb_mode_S,nb_mode_S) )

KSS_static_S1 = scipy.dot(KSS[SolvedDofS,:][:,SolvedDofS],scipy.sparse.coo_matrix(S).T)
staticT_KSS_static_11 = scipy.array(S*KSS_static_S1.todense())[0][0]

Knn = scipy.sparse.construct.bmat( [[K_diag_nn,None],[None,staticT_KSS_static_11]] )

Mnn = scipy.sparse.construct.bmat( [[M_diag_nn,None],[None,1.0]] )



K=scipy.sparse.construct.bmat( [ [fluid_damping*K_diag_pp,K_pk,None],[K_pk.T,fluid_damping*K_kk,None],[None,-Cnk,Knn] ] )

M=scipy.sparse.construct.bmat( [ [M_diag_pp,M_pk,None],[M_pk.T,M_kk,Cnk.T],[None,None,Mnn] ] )

F  = scipy.zeros((nb_mode_F+nb_mode_A))
Fn = PSn.T*FS[SolvedDofS]
F  = scipy.append(F,Fn)
F  = scipy.sparse.csc_matrix(F).T

##############################################################

if rank==0:
    print ("nb. fluid modes = ",nb_mode_F)
    print ("nb. enriched nodes = ",len(SolvedDofA))
    print ("nb. structure modes = ",nb_mode_S)

##############################################################
# FRF computation of the FSI problem
##############################################################
toc0=time.clock()
time_before_frf=time.ctime()
if rank==0:
    print ("fixed time before FRF loop:",toc0-tic0)


Flag_frf_analysis=1
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    if rank==0:
        print ("time at the beginning of the FRF:",time.ctime())

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):
        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq
        if rank==0:
            print ("proc number",rank,"frequency=",freq)

        #sol = scipy.sparse.linalg.spsolve(K-(omega*omega)*M+omega*D*1j, F)
        #sol = scipy.linalg.solve(scipy.array((K-(omega*omega)*M+omega*D*1j).todense()),scipy.array(F.todense()))
        sol = mumps.spsolve(  scipy.sparse.csc_matrix(K-(omega*omega)*M+omega*D*1j,dtype='c16')  , scipy.array(F.todense() , dtype='c16'), comm=mycomm )

        alpha_p    = scipy.sparse.csc_matrix(sol[list(range(nb_mode_F))])
        alpha_k    = scipy.sparse.csc_matrix(sol[list(range(nb_mode_F,nb_mode_F+nb_mode_A,1))])
        P_A        = Psi_Ak*alpha_k
        press      = eigen_vectors_F*alpha_p+Psi_Fk*alpha_k
        enrichment = scipy.zeros((fluid_nnodes),dtype=complex)
        enrichment[SolvedDofA]=P_A
        CorrectedPressure=press
        CorrectedPressure[SolvedDofA]=(CorrectedPressure[SolvedDofA].T+scipy.array(enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)).T
        frf.append(silex_lib_xfem_acou_tet4.computecomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
        #frf.append(silex_lib_xfem_acou_tet4.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press.todense(),enrichment,LevelSet,LevelSetTangent))

        if rank==0:
            Q=scipy.zeros((struc_ndof),dtype=float)
            tmp=sol[list(range(nb_mode_F+nb_mode_A,nb_mode_F+nb_mode_A+PSn.shape[1],1))].real
            Q[SolvedDofS]=PSn*tmp
            disp=scipy.zeros((struc_nnodes,3))
            disp[range(struc_nnodes),0]=Q[list(range(0,struc_ndof,6))]
            disp[range(struc_nnodes),1]=Q[list(range(1,struc_ndof,6))]
            disp[range(struc_nnodes),2]=Q[list(range(2,struc_ndof,6))]
            disp_save.append(disp)
            press_save.append(CorrectedPressure.real)

    if rank==0:
        print ("Proc. ",rank," / time at the end of the FRF:",time.ctime())
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements,4,[[press_save,'nodal',1,'pressure']])
        silex_lib_gmsh.WriteResults2(results_file+str(rank)+'_results_struct_frf',struc_nodes,struc_elements,2,[[disp_save,'nodal',3,'displacement']])
    
    # Save the FRF problem
    Allfrequencies=scipy.zeros(nb_freq_step)
    Allfrf=scipy.zeros(nb_freq_step)
    k=0
    if rank==0:
        for i in range(nproc):
            data = comm.recv(source=i, tag=11)
            for j in range(len(data[0])):
                Allfrequencies[k]=data[0][j]
                Allfrf[k]=data[1][j]
                k=k+1

        Allfrequencies, Allfrf = zip(*sorted(zip(Allfrequencies, Allfrf)))
        Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf))]
        f=open(results_file+'_results.frf','wb')
        pickle.dump(Allfrfsave, f)
        f.close()

        print("rank before Gram-Schmidt =",numpy.linalg.matrix_rank(Psi))
        print("rank after  Gram-Schmidt =",numpy.linalg.matrix_rank(Psi_ortho))
        print("rank of structure basis  =",numpy.linalg.matrix_rank(scipy.array(PSn.todense())))

        #print ("structure eigen frequencies (fixed): ",freq_eigv_S)
        #print ("-------")
        #print ("fluid eigen frequencies : ",freq_eigv_F)

        if rank==0:
            print("last freq / fluid / ",freq_fluid[nb_mode_F-1])
            print("last freq / struc / ",freq_struc_fixed[nb_mode_S-1])

            if nb_mode_S_enc_k!=0:
                print("last freq / interface, fixed modes / ",freq_struc_fixed[nb_mode_S_enc_k-1])
            else:
                print("last freq / interface, fixed modes / None")
                
            if nb_mode_S_fre_k!=0:
                print("last freq / interface, free  modes / ",freq_struc_free[nb_mode_S_fre_k-1])
            else:
                print("last freq / interface, free  modes / None")

        print("Size of the reduced problem = ",K.shape)

        print("Real time at the beginning = ",time_init)
        print("Real time before FRF       = ",time_before_frf)
        print("Real time at the end       = ",time.ctime())
        print("Total time = ",time.clock()-tic00)

