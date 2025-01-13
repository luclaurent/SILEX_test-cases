import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle

#import mumps

import sys
sys.path.append('../librairies')
import silex_lib_gmsh
import silex_acou_lib_tri3

from mpi4py import MPI
comm = MPI.COMM_WORLD

nproc=comm.Get_size()
rank = comm.Get_rank()

# mpirun -np 2 python Main_xfem.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
#mesh_file='results/xfem1';results_file='results/xfem-test1-no-edge'
#mesh_file='results/xfem2';results_file='results/xfem-test2-no-edge'
#mesh_file='results/xfem3';results_file='results/xfem-test3-no-edge'
#mesh_file='results/xfem4';results_file='results/xfem-test4-no-edge'
#mesh_file='results/xfem5';results_file='results/xfem-test5-no-edge'
mesh_file='results/xfem6';results_file='results/xfem-test6-no-edge'

##############################################################
# Material, Boundary conditions
##############################################################

# air
celerity=340.0
rho=1.2

freq_ini     = 100.5
freq_end     = 300.0
nb_freq_step_per_proc=50 # 50 pour 8 proc.

nproc=comm.Get_size()
rank = comm.Get_rank()

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=1

flag_edge_enrichment=0
#flag_edge_enrichment=1

freq_comparaison = 210.0


##############################################################
# Load fluid mesh
##############################################################

tic = time.clock()

fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_fluid.msh',2)
fluid_elements = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',2,1)

fluid_nnodes   = fluid_nodes.shape[0]
fluid_nelem    = fluid_elements.shape[0]
fluid_ndof     = fluid_nnodes

fluid_elements_boun = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',1,2)
IdnodeS2=scipy.unique(fluid_elements_boun)

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,2)
    silex_lib_gmsh.WriteResults(results_file+'_fluid_boundary',fluid_nodes,fluid_elements_boun,1)
print 'nnodes for fluid=',fluid_nnodes
print 'nelem for fluid=',fluid_nelem



##################################################################
# compute level set
##################################################################

tic = time.clock()

LevelSet=fluid_nodes[:,0]-0.6
LevelSetTangent=fluid_nodes[:,1]-0.65

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,2,[[[LevelSetTangent],'nodal',1,'Tangent level set']])
    silex_lib_gmsh.WriteResults(results_file+'_level_set',fluid_nodes,fluid_elements,2,[[[LevelSet],'nodal',1,'Level set']])

toc = time.clock()
print "time to compute level set:",toc-tic


##################################################################
# Get enriched nodes and elements
##################################################################
tic = time.clock()

struc_nodes=scipy.array([[0.6,0.0],[0.6,0.3],[0.6,0.65]])
struc_elements=scipy.array([[1,2],[2,3]])
struc_boun=scipy.array([3])

silex_lib_gmsh.WriteResults(results_file+'_struc_mesh',struc_nodes,struc_elements,1)

EnrichedElements,NbEnrichedElements=silex_acou_lib_tri3.getenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements)
EnrichedElements=scipy.unique(EnrichedElements[range(NbEnrichedElements)])-1

toc = time.clock()
print "time to find surface enriched elements:",toc-tic

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],2)

tic = time.clock()

EdgeEnrichedElements,nbenrelts = silex_acou_lib_tri3.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes,fluid_elements)
EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[range(nbenrelts)])-1

toc = time.clock()
print "time to find edge enriched elements:",toc-tic

if (flag_write_gmsh_results==1) and (rank==0):
    silex_lib_gmsh.WriteResults(results_file+'_edge_enriched_elements',fluid_nodes,fluid_elements[EdgeEnrichedElements],2)

##############################################################
# Compute Standard Fluid Matrices
##############################################################
tic = time.clock()

IIf,JJf,Vffk,Vffm=silex_acou_lib_tri3.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

SolvedDofF=scipy.setdiff1d(range(fluid_ndof),IdnodeS2-1)
#SolvedDofF=range(fluid_ndof)

toc = time.clock()
print "time to compute fluid matrices:",toc-tic

##################################################################
# Compute enrichment: Heaviside + Edge
##################################################################
tic = time.clock()

#HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

#Enrichednodes = scipy.unique(fluid_elements[HeavisideEnrichedElements])
#Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

#print xvibacoufo.getpositivenegativeelts.__doc__

NegativeLSelements,PositiveLSelements,NegativeLStgtElements,PositiveLStgtElements,nbNegLS,nbPosLS,nbNegLSt,nbPosLSt=silex_acou_lib_tri3.getpositivenegativeelts(fluid_elements,LevelSet,LevelSetTangent)

NegativeLSelements=NegativeLSelements[range(nbNegLS)]
PositiveLSelements=PositiveLSelements[range(nbPosLS)]
NegativeLStgtElements=NegativeLStgtElements[range(nbNegLSt)]
PositiveLStgtElements=PositiveLStgtElements[range(nbPosLSt)]

EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_acou_lib_tri3.getenrichedelementsfromlevelset(fluid_elements,LevelSetTangent)
EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[range(nbEdgeEnrichedElementsInAllMesh)])

IdElementTip=silex_acou_lib_tri3.getelementcontainingpoint(fluid_elements,fluid_nodes,[0.6,0.65])

#if (flag_write_gmsh_results==1) and (rank==0):
#    silex_lib_gmsh.WriteResults(results_file+'_NegativeLSelements',fluid_nodes,fluid_elements[NegativeLSelements],2)
#    silex_lib_gmsh.WriteResults(results_file+'_PositiveLSelements',fluid_nodes,fluid_elements[PositiveLSelements],2)
#    silex_lib_gmsh.WriteResults(results_file+'_NegativeLStgtElements',fluid_nodes,fluid_elements[NegativeLStgtElements],2)
#    silex_lib_gmsh.WriteResults(results_file+'_PositiveLStgtElements',fluid_nodes,fluid_elements[PositiveLStgtElements],2)




IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_acou_lib_tri3.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,LevelSetTangent,celerity,rho,flag_edge_enrichment)

KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

toc = time.clock()
print 'time to compute Heaviside enrichment:',toc-tic

#Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([HeavisideEnrichedElements,EdgeEnrichedElements]))])
#Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements,EdgeEnrichedElementsInAllMesh]))])
#Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([EnrichedElements,PositiveLStgtElements]))])
#Enrichednodes = scipy.unique(fluid_elements[scipy.hstack(([NegativeLStgtElements]))])
Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])
#Enrichednodes = scipy.unique(fluid_elements)
SolvedDofA=Enrichednodes-1

silex_lib_gmsh.WriteResults(results_file+'_EnrichedElements',fluid_nodes,fluid_elements[scipy.hstack(([EnrichedElements]))],2)

#################################################################
# Construct the whole system
##################################################################

K=scipy.sparse.construct.bmat( [[KFF[scipy.ix_(SolvedDofF,SolvedDofF)],KAF[scipy.ix_(SolvedDofF,SolvedDofA)]],
                                [KAF[scipy.ix_(SolvedDofA,SolvedDofF)],KAA[scipy.ix_(SolvedDofA,SolvedDofA)]]
                                ] )



M=scipy.sparse.construct.bmat( [[MFF[scipy.ix_(SolvedDofF,SolvedDofF)],MAF[scipy.ix_(SolvedDofF,SolvedDofA)]],
                                [MAF[scipy.ix_(SolvedDofA,SolvedDofF)],MAA[scipy.ix_(SolvedDofA,SolvedDofA)]]
                                ] )

##############################################################
# FRF computation of the FSI problem
##############################################################

Flag_frf_analysis=1
FF = scipy.zeros(fluid_ndof)
frequencies=[]
frf=[]

if (Flag_frf_analysis==1):
    print "time at the beginning of the FRF:",time.ctime()

    press_save=[]
    disp_save=[]

    for i in range(nb_freq_step_per_proc):

        freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
        frequencies.append(freq)
        omega=2*scipy.pi*freq
        print "proc number",rank,"frequency=",freq

        FF[scipy.ix_(SolvedDofF)]=-(KFF[scipy.ix_(SolvedDofF,IdnodeS2-1)]-(omega*omega)*MFF[scipy.ix_(SolvedDofF,IdnodeS2-1)])*(scipy.zeros((len(IdnodeS2)))+1.0)
        FA = scipy.zeros(fluid_ndof)
        F  = FF[scipy.ix_(SolvedDofF)]
        F  = scipy.append(F,FA[scipy.ix_(SolvedDofA)])
        F  = scipy.sparse.csc_matrix(F)

        sol = scipy.sparse.linalg.spsolve(K-(omega*omega)*M, F)

        press = scipy.zeros(fluid_ndof)
        press[scipy.ix_(SolvedDofF)]=sol[range(len(SolvedDofF))]
        enrichment=scipy.zeros(fluid_nnodes)
        enrichment[scipy.ix_(SolvedDofA)]=sol[range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA))]
        CorrectedPressure=press
        CorrectedPressure[scipy.ix_(SolvedDofA)]=CorrectedPressure[scipy.ix_(SolvedDofA)]+enrichment[scipy.ix_(SolvedDofA)]*scipy.sign(LevelSet[scipy.ix_(SolvedDofA)])
        #frf.append(silex_acou_lib_tri3.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
        frf.append(silex_acou_lib_tri3.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,press+0j,enrichment+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
        #frf[i]=xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent)
        press_save.append(CorrectedPressure)
        
        if freq==freq_comparaison:
            PressureAtTip,ThetaAtTip=silex_acou_lib_tri3.computexfempressureattip(fluid_elements,fluid_nodes,press+0j,enrichment+0j,[0.6,0.65])

            AllPressTip=[PressureAtTip,ThetaAtTip,freq]
            f=open(results_file+'_pressTip.frf','w')
            pickle.dump(AllPressTip, f)
            f.close()


    print "time at the end of the FRF:",time.ctime()
    frfsave=[frequencies,frf]
    comm.send(frfsave, dest=0, tag=11)
    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+str(rank)+'_results_fluid_frf',fluid_nodes,fluid_elements,2,[[press_save,'nodal',1,'pressure']])

    # save the FRF problem
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
        f=open(results_file+'_results.frf','w')
        pickle.dump(Allfrfsave, f)
        f.close()

    


