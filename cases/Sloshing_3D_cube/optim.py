import Main_parameters_tet4 as solver
from scipy import optimize

obj = solver.solve()

val_bnd = 0.21
bounds = [(-0.21, 0.21), (-0.21, 0.21)]

results = dict()
results['shgo'] = optimize.shgo(obj.run, bounds)