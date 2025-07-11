from pathlib import Path
from meshRW import msh, msh2, vtk
import numpy as np
from SILEXlib import MeshField as lib
from SILEXlib import silex_lib_xfem_acou_tet10 as libF_tet10_xfem


# mesh file
mesh_file = Path(__file__).parent / 'debug_tet10.msh'
results_file = Path(__file__).parent / 'debug_tet10_'

# LS position
pos = 0.25

# read mesh file
mesh = msh.mshReader(mesh_file)
nodes = mesh.getNodes()
elements = mesh.getElements(type='TET10')

# convert TET10 to TET4
# elementsTET4fromTET10 = elements[:,0:4]
elementsTET4fromTET10 = libF_tet10_xfem.tet10totet4(nodes, elements)

# create LS
LS = nodes[:,0] - pos
LST = np.zeros_like(LS)-1.0

# create pressures fields
uncorrected_pressure = np.zeros(nodes.shape[0])
enrich_pressure = np.zeros(nodes.shape[0])
# uncorrected_pressure[:] = 1.0
# enrich_pressure[:] = 1.0

k = 1
lon = 1.1
enrich_pressure[:] = np.sin(nodes[:,0]*k*np.pi/ lon)    
uncorrected_pressure[:] = np.sin(nodes[:,0]*k*np.pi/ lon)   



# run post-processing utility
objMesh = lib.MeshField(nodes=nodes, 
                        elems=elementsTET4fromTET10, 
                        leveset=LS, 
                        levelsetTg=LST)
objMesh.addField(uncorrectedField=uncorrected_pressure,
                enrichmentField=enrich_pressure)
datameshfield = objMesh.getData()

# write results file
# prepare fields
dataW = []
dataW.append({'data':datameshfield['LS'],'type':'nodal', 'name':'levelset'})
dataW.append({'data':datameshfield['LST'],'type':'nodal','name':'tangent levelset'})
# for it in range(len(press)):
#     dataW.append({'data':datameshfield['fields'][:,it],'type':'nodal','name':'field '+str(it)+' (levelset)'})
dataW.append({
    'name': 'press',
    'type': 'nodal',
    'data': datameshfield['fields']
})
dataW.append({
    'name': 'uncorrected press',
    'type': 'nodal',
    'data': datameshfield['uncorrected']
})
dataW.append({
    'name': 'correction press',
    'type': 'nodal',
    'data': datameshfield['correction']
})
#            # export mesh
msh2.mshWriter(
    filename= results_file.as_posix() +'_results_fluid_frf_meshfield3D.msh',
    nodes=datameshfield['nodes'],
    elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
            {'type':'PRI6','connectivity':datameshfield['PRI6']}],
    fields=dataW,
    append=True
    )

# vtk.vtkWriter(
#     filename= results_file.as_posix() +'_results_fluid_frf_meshfield3D.vtk',
#     nodes=datameshfield['nodes'],
#     elements=[{'type':'TET4','connectivity':datameshfield['TET4']},
#             {'type':'PRI6','connectivity':datameshfield['PRI6']}],
#     fields=dataW,
#     append=True
#     )



