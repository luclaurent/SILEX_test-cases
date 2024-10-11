import string
import time
import scipy
import scipy.sparse
import scipy.sparse.linalg

import pylab as pl
import pickle

#import mumps

import sys
sys.path.append('../../librairies')
import silex_lib_gmsh
import silex_lib_tri3_acou

from mpi4py import MPI
comm = MPI.COMM_WORLD

nproc=comm.Get_size()
rank = comm.Get_rank()

import mumps

from mpi4py import MPI
comm = MPI.COMM_WORLD
nproc=comm.Get_size()
rank = comm.Get_rank()

class comm_mumps_one_proc:
    rank = 0
    def py2f(self):
        return 0
mycomm=comm_mumps_one_proc()

###########################################################
# To run it in parallel for several frequencies:
# export OPENBLAS_NUM_THREADS=1
# mpirun -np 20 python3 Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
#
# To run it in sequentiel frequency per frequency with openblas in parrallel:
# export OPENBLAS_NUM_THREADS=10
# python3 Main_toto.py

# mpirun -np 2 python Main_xfem_2.py

##############################################################
##############################################################
#              S T A R T   M A I N   P R O B L E M
##############################################################
##############################################################

# parallepipedic cavity with plane structure
mesh_file='xfem3'
results_fileInit='xfem_3_'

##############################################################
# Material, Boundary conditions
##############################################################

# air
celerity=340.0
rho=1.2

freq_ini     = 100.0
freq_end     = 300.0
nb_freq_step_per_proc=50 # 50 pour 8 proc.

nproc=comm.Get_size()
rank = comm.Get_rank()

nb_freq_step = nb_freq_step_per_proc*nproc
deltafreq=(freq_end-freq_ini)/(nb_freq_step-1)

flag_write_gmsh_results=1

flag_edge_enrichment=0
#flag_edge_enrichment=1

#freq_comparaison = 210.0

for itP in [0,1,5,10,50,100,500]:

    x_pos_struc=620000+itP
    h_struc=1.0

    file_extension=str(x_pos_struc)[0:6]
    results_file=results_fileInit+file_extension

    x_pos_struc=x_pos_struc/1000000

    ##############################################################
    # Load fluid mesh
    ##############################################################

    tic = time.clock()

    fluid_nodes    = silex_lib_gmsh.ReadGmshNodes(mesh_file+'_fluid.msh',2)
    fluid_elements,Idnodes = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',2,1)
    fluid_elements5,Idnodes5 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',2,5)

    fluid_nnodes   = fluid_nodes.shape[0]
    fluid_nelem    = fluid_elements.shape[0]
    fluid_ndof     = fluid_nnodes

    fluid_elements_boun,IdnodeS2 = silex_lib_gmsh.ReadGmshElements(mesh_file+'_fluid.msh',1,2)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_fluid_mesh',fluid_nodes,fluid_elements,2)
        silex_lib_gmsh.WriteResults(results_file+'_control_volume_fluid_mesh',fluid_nodes,fluid_elements5,2)
        silex_lib_gmsh.WriteResults(results_file+'_fluid_boundary',fluid_nodes,fluid_elements_boun,1)
    print("nnodes for fluid=",fluid_nnodes)
    print("nelem for fluid=",fluid_nelem)



    ##################################################################
    # compute level set
    ##################################################################

    tic = time.clock()

    LevelSet=fluid_nodes[:,0]-x_pos_struc
    LevelSetTangent=fluid_nodes[:,1]-h_struc

    # level set gradient with respect to parameters
    LevelSet_gradient=-scipy.ones(fluid_nnodes)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_tangent_level_set',fluid_nodes,fluid_elements,2,[[LevelSetTangent,'nodal',1,'Tangent level set']])
        silex_lib_gmsh.WriteResults(results_file+'_level_set',fluid_nodes,fluid_elements,2,[[LevelSet,'nodal',1,'Level set']])

    toc = time.clock()
    print("time to compute level set:",toc-tic)


    ##################################################################
    # Get enriched nodes and elements
    ##################################################################
    tic = time.clock()

    struc_nodes=scipy.array([[x_pos_struc,0.0],[x_pos_struc,h_struc/2.0],[x_pos_struc,h_struc]])
    struc_elements=scipy.array([[1,2],[2,3]])
    struc_boun=scipy.array([3])

    silex_lib_gmsh.WriteResults(results_file+'_struc_mesh',struc_nodes,struc_elements,1)

    EnrichedElements,NbEnrichedElements=silex_lib_tri3_acou.getenrichedelements(struc_nodes,struc_elements,fluid_nodes,fluid_elements)
    EnrichedElements=scipy.unique(EnrichedElements[list(range(NbEnrichedElements))])-1

    toc = time.clock()
    print("time to find surface enriched elements:",toc-tic)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_enriched_elements',fluid_nodes,fluid_elements[EnrichedElements],2)

    tic = time.clock()

    EdgeEnrichedElements,nbenrelts = silex_lib_tri3_acou.getedgeenrichedelements(struc_nodes,struc_boun,fluid_nodes,fluid_elements)
    EdgeEnrichedElements=scipy.unique(EdgeEnrichedElements[list(range(nbenrelts))])-1

    toc = time.clock()
    print("time to find edge enriched elements:",toc-tic)

    if (flag_write_gmsh_results==1) and (rank==0):
        silex_lib_gmsh.WriteResults(results_file+'_edge_enriched_elements',fluid_nodes,fluid_elements[EdgeEnrichedElements],2)

    ##############################################################
    # Compute Standard Fluid Matrices
    ##############################################################
    tic = time.clock()

    IIf,JJf,Vffk,Vffm=silex_lib_tri3_acou.globalacousticmatrices(fluid_elements,fluid_nodes,celerity,rho)

    KFF = scipy.sparse.csc_matrix( (Vffk,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    MFF = scipy.sparse.csc_matrix( (Vffm,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    SolvedDofF=scipy.setdiff1d(list(range(fluid_ndof)),IdnodeS2-1)
    #SolvedDofF=range(fluid_ndof)

    toc = time.clock()
    print("time to compute fluid matrices:",toc-tic)

    ##################################################################
    # Compute enrichment: Heaviside + Edge
    ##################################################################
    tic = time.clock()

    #HeavisideEnrichedElements=scipy.setdiff1d(EnrichedElements,EdgeEnrichedElements)

    #Enrichednodes = scipy.unique(fluid_elements[HeavisideEnrichedElements])
    #Enrichednodes = scipy.unique(fluid_elements[EnrichedElements])

    #print xvibacoufo.getpositivenegativeelts.__doc__

    NegativeLSelements,PositiveLSelements,NegativeLStgtElements,PositiveLStgtElements,nbNegLS,nbPosLS,nbNegLSt,nbPosLSt=silex_lib_tri3_acou.getpositivenegativeelts(fluid_elements,LevelSet,LevelSetTangent)

    NegativeLSelements=NegativeLSelements[list(range(nbNegLS))]
    PositiveLSelements=PositiveLSelements[list(range(nbPosLS))]
    NegativeLStgtElements=NegativeLStgtElements[list(range(nbNegLSt))]
    PositiveLStgtElements=PositiveLStgtElements[list(range(nbPosLSt))]

    EdgeEnrichedElementsInAllMesh,nbEdgeEnrichedElementsInAllMesh=silex_lib_tri3_acou.getenrichedelementsfromlevelset(fluid_elements,LevelSetTangent)
    EdgeEnrichedElementsInAllMesh=scipy.unique(EdgeEnrichedElementsInAllMesh[list(range(nbEdgeEnrichedElementsInAllMesh))])

    IdElementTip=silex_lib_tri3_acou.getelementcontainingpoint(fluid_elements,fluid_nodes,[0.6,0.65])

    #if (flag_write_gmsh_results==1) and (rank==0):
    #    silex_lib_gmsh.WriteResults(results_file+'_NegativeLSelements',fluid_nodes,fluid_elements[NegativeLSelements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_PositiveLSelements',fluid_nodes,fluid_elements[PositiveLSelements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_NegativeLStgtElements',fluid_nodes,fluid_elements[NegativeLStgtElements],2)
    #    silex_lib_gmsh.WriteResults(results_file+'_PositiveLStgtElements',fluid_nodes,fluid_elements[PositiveLStgtElements],2)




    IIaa,JJaa,IIaf,JJaf,Vaak,Vaam,Vafk,Vafm=silex_lib_tri3_acou.globalxfemacousticmatrices(fluid_elements,fluid_nodes,LevelSet,LevelSetTangent,celerity,rho,flag_edge_enrichment)

    KAA = scipy.sparse.csc_matrix( (Vaak,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    MAA = scipy.sparse.csc_matrix( (Vaam,(IIaa,JJaa)), shape=(fluid_ndof,fluid_ndof) )
    KAF = scipy.sparse.csc_matrix( (Vafk,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )
    MAF = scipy.sparse.csc_matrix( (Vafm,(IIaf,JJaf)), shape=(fluid_ndof,fluid_ndof) )

    toc = time.clock()
    print("time to compute Heaviside enrichment:",toc-tic)

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

    K=scipy.sparse.construct.bmat( [[KFF[SolvedDofF,:][:,SolvedDofF],KAF[SolvedDofF,:][:,SolvedDofA]],
                                    [KAF[SolvedDofA,:][:,SolvedDofF],KAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )



    M=scipy.sparse.construct.bmat( [[MFF[SolvedDofF,:][:,SolvedDofF],MAF[SolvedDofF,:][:,SolvedDofA]],
                                    [MAF[SolvedDofA,:][:,SolvedDofF],MAA[SolvedDofA,:][:,SolvedDofA]]
                                    ] )

    #################################################################
    # Compute gradients with respect to parameters
    ##################################################################

    print(silex_lib_tri3_acou.globalacousticgradientmatrices.__doc__)

    IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements[EnrichedElements],fluid_nodes,celerity,rho,LevelSet_gradient,LevelSet)
    #IIf,JJf,Vfak_gradient,Vfam_gradient=silex_lib_tri3_acou.globalacousticgradientmatrices(fluid_elements,fluid_nodes,celerity,rho,LevelSet_gradient,LevelSet)
    dMFA_dtheta = scipy.sparse.csc_matrix( (Vfam_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )
    dKFA_dtheta = scipy.sparse.csc_matrix( (Vfak_gradient,(IIf,JJf)), shape=(fluid_ndof,fluid_ndof) )

    f=open(results_file+'_KAF_MAF.pck','wb')
    pickle.dump([dKFA_dtheta[SolvedDofF,:][:,SolvedDofA],dMFA_dtheta[SolvedDofF,:][:,SolvedDofA],KAF[SolvedDofF,:][:,SolvedDofA],MAF[SolvedDofF,:][:,SolvedDofA]], f)
    f.close()

    dK=scipy.sparse.construct.bmat( [[None,dKFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [dKFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )

    dM=scipy.sparse.construct.bmat( [[None,dMFA_dtheta[SolvedDofF,:][:,SolvedDofA]],
                                     [dMFA_dtheta[SolvedDofA,:][:,SolvedDofF],None]
                                    ] )


    ##############################################################
    # FRF computation of the FSI problem
    ##############################################################

    Flag_frf_analysis=1
    FF = scipy.zeros(fluid_ndof)
    frequencies=[]
    frf=[]
    frfgradient=[]

    if (Flag_frf_analysis==1):
        print("time at the beginning of the FRF:",time.ctime())

        press_save=[]
        dpress_save=[]
        disp_save=[]

        for i in range(nb_freq_step_per_proc):

            freq = freq_ini+i*nproc*deltafreq+rank*deltafreq
            frequencies.append(freq)
            omega=2*scipy.pi*freq
            print("proc number",rank,"frequency=",freq)

            FF[SolvedDofF]=-(KFF[SolvedDofF,:][:,IdnodeS2-1]-(omega**2)*MFF[SolvedDofF,:][:,IdnodeS2-1])*(scipy.ones((len(IdnodeS2))))
            FA = scipy.zeros(fluid_ndof)
            F  = FF[SolvedDofF]
            F  = scipy.append(F,FA[SolvedDofA])
            #F  = scipy.sparse.csc_matrix(F)

            sol = mumps.spsolve(scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float'), F, comm=mycomm )
            
            tmp=-(dK-(omega**2)*dM)*sol
            Dsol_Dtheta = mumps.spsolve(  scipy.sparse.coo_matrix(K-(omega**2)*M,dtype='float')  , tmp , comm=mycomm )


            press = scipy.zeros(fluid_ndof)
            press[IdnodeS2-1] = scipy.ones(len(IdnodeS2))
            press[SolvedDofF]=sol[list(range(len(SolvedDofF)))]
            enrichment=scipy.zeros(fluid_nnodes)
            enrichment[SolvedDofA]=sol[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            CorrectedPressure=press
            CorrectedPressure[SolvedDofA]=CorrectedPressure[SolvedDofA]+enrichment[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA])
            #frf.append(silex_acou_lib_tri3.computequadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure))
            frf.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressure(fluid_elements5,fluid_nodes,press+0j,enrichment+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            #frf[i]=xvibacoufo.computexfemcomplexquadratiquepressure(fluid_elements,fluid_nodes,CorrectedPressure+0j,0.0*enrichment+0j,LevelSet,LevelSetTangent)
            press_save.append(CorrectedPressure)

            Dpress_Dtheta = scipy.zeros(fluid_ndof,dtype=float)
            Dpress_Dtheta[SolvedDofF] = Dsol_Dtheta[list(range(len(SolvedDofF)))]
            Denrichment_Dtheta = scipy.zeros(fluid_ndof,dtype=float)
            Denrichment_Dtheta[SolvedDofA]= Dsol_Dtheta[list(range(len(SolvedDofF),len(SolvedDofF)+len(SolvedDofA)))]
            #DCorrectedPressure_Dtheta=scipy.array(Dpress_Dtheta)
            #DCorrectedPressure_Dtheta[SolvedDofA]=DCorrectedPressure_Dtheta[SolvedDofA].T+scipy.array(Denrichment_Dtheta[SolvedDofA]*scipy.sign(LevelSet[SolvedDofA]).T)
            frfgradient.append(silex_lib_tri3_acou.computexfemcomplexquadratiquepressuregradient(fluid_elements5,fluid_nodes,press+0j,enrichment+0j,Dpress_Dtheta+0j,Denrichment_Dtheta+0j,LevelSet,LevelSetTangent,flag_edge_enrichment))
            dpress_save.append(Dpress_Dtheta)
        

        print("time at the end of the FRF:",time.ctime())
        frfsave=[frequencies,frf]
        comm.send(frfsave, dest=0, tag=11)
        if (flag_write_gmsh_results==1) and (rank==0):
            silex_lib_gmsh.WriteResults2(results_file+'_results_fluid_frf',fluid_nodes,fluid_elements,2,[[press_save,'nodal',1,'pressure'],[dpress_save,'nodal',1,'pressure gradient']])

        # save the FRF problem
        Allfrequencies=scipy.zeros(nb_freq_step)
        Allfrf=scipy.zeros(nb_freq_step)
        Allfrfgradient=scipy.zeros(nb_freq_step)
        k=0
        if rank==0:
            for i in range(nproc):
                data = comm.recv(source=i, tag=11)
                for j in range(len(data[0])):
                    Allfrequencies[k]=data[0][j]
                    Allfrf[k]=data[1][j]
                    Allfrfgradient[k]=data[1][j]
                    k=k+1

            Allfrequencies, Allfrf,Allfrfgradient = zip(*sorted(zip(Allfrequencies, Allfrf,Allfrfgradient)))
            Allfrfsave=[scipy.array(list(Allfrequencies)),scipy.array(list(Allfrf)),scipy.array(list(Allfrfgradient)),press_save,dpress_save]
            f=open(results_file+'_results.frf','wb')
            pickle.dump(Allfrfsave, f)
            f.close()


