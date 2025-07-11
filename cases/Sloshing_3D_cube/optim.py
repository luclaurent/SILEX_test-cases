import Main_parameters_tet4 as solver
from scipy import optimize

obj = solver.solve()

bounds = [(-0.29, 0.29), (-0.29, 0.29)]

results = dict()
results['shgo'] = optimize.shgo(obj.run, bounds)