#!/usr/bin/env python
# coding: utf-8

# ### run CBO on 6 parameters case

# In[2]:




# In[ ]:


import numpy as np
# build mechanical function
from Main_parameters_tet4 import solve
class solver:
    def __init__(self, nb_para=2):
        self.X = np.zeros((0,nb_para))
        self.obj = np.zeros(0)
        self.cons = np.zeros(0)
        self.model = solve()
        self.fmax = 19.0
    def formatX(self,x):
        if len(x.shape)==1:
            return x
        elif len(x.shape)==2 and x.shape[0]==1:
            return x[0]
    def run(self,x):
        x_formatted = self.formatX(x)
        val_p,_, val_f,_,_,_ = self.model.run_bis(x_formatted)
        self.X = np.vstack((self.X, x_formatted))
        self.obj = np.hstack((self.obj, val_p/(1000*9.81)))
        cons_val = self.fmax-val_f
        self.cons = np.hstack((self.cons, cons_val))
        return val_p, cons_val
    def fun_obj(self,x):
        x_formatted = self.formatX(x)
        IXraw = np.where((self.X==x_formatted).all(axis=1))[0]
        if len(IXraw) > 0:
            return self.obj[IXraw[0]]
        val_p, _ = self.run(x_formatted)  
        return val_p
    def fun_cons(self,x):
        x_formatted = self.formatX(x)
        IXraw = np.where((self.X==x_formatted).all(axis=1))[0]
        if len(IXraw) > 0:
            return self.cons[IXraw[0]]
        _, val_c = self.run(x_formatted)  
        return val_c
    
mechanical_model = solver(nb_para=6)
    


# In[4]:


bounds_6D = [[-0.2, -0.2,-0.2, -0.2,-0.1, -0.1], 
          [0.2, 0.2,0.2, 0.2,0.1, 0.1]]





# ### run complete CBO study

# In[ ]:

import bo_lib
import torch
from loguru import logger
nb_runs = 5
nb_it_cbo = 50
nb_samples = [5,10,15]
all_data = list()
data_for_ns = list()
for itns,ns in enumerate(nb_samples):
    data_for_ns.append(list())
    for _ in range(nb_runs):
        # sampling
        samples_6D = bo_lib.lhs_distrib(bounds_6D, nbs=ns)
        Zsamples_6D = np.zeros((0,))
        Csamples_6D = np.zeros((0,))
        for p in samples_6D:
            Ztmp=mechanical_model.fun_obj(p)
            Ctmp=mechanical_model.fun_cons(p)
            Zsamples_6D = np.hstack((Zsamples_6D, Ztmp))
            Csamples_6D = np.hstack((Csamples_6D, Ctmp))
        Zsamples_6D = torch.tensor(Zsamples_6D[np.newaxis].T, dtype=torch.float64)
        Csamples_6D = torch.tensor(Csamples_6D[np.newaxis].T, dtype=torch.float64)


        dataStdizeX_6D = bo_lib.stdize(samples_6D)
        dataStdizeZ_6D = bo_lib.stdize(Zsamples_6D)
        dataStdizeC_6D = bo_lib.stdize(Csamples_6D)
        Xstdize_6D = dataStdizeX_6D.output_data
        Zstdize_6D = dataStdizeZ_6D.output_data
        Cstdize_6D = dataStdizeC_6D.output_data
        Cstdize_6D = [Cstdize_6D]
        # build multi independent GP
        multiZC_6D = torch.hstack((Zstdize_6D, *Cstdize_6D))
        gp_train_multi_6D = bo_lib.init_surrogate_model(Xstdize_6D, multiZC_6D)     

        # check interpolation
        multi_GP_init_at_samples_6D = gp_train_multi_6D(Xstdize_6D).mean.detach()
        for itGP,gpI in enumerate(multi_GP_init_at_samples_6D):
            logger.info(f'Index {itGP} maxi interpolation error {(multiZC_6D[:,itGP].flatten()-gpI).abs().max()}')
        # run BO
        constraints_bounds_6D = [(0,torch.inf)]
        constraints_bounds_stdize_6D = [dataStdizeC_6D.stdize(torch.tensor(bnd, dtype=torch.float64)) for bnd in constraints_bounds_6D]
        normalization_data_6D={
            'X': dataStdizeX_6D,
            'Z': dataStdizeZ_6D,
            'C': dataStdizeC_6D
        }
        # run BO
        (
            CBO_X_update_6D, 
            CBO_ZC_update_6D, 
            CBO_gp_obj_final_6D,     
            CBO_GP_obj_at_grid_points_6D, 
            CBO_GP_cons_at_grid_points_6D,
            CBO_acq_at_grid_points_6D,
            acq_values_final_6D,
            best_candidates,
            X_best_candidates
            ) = bo_lib.execute_CBO(gp_train_multi_6D,
                            fun_obj=mechanical_model.fun_obj,
                            fun_cons=[mechanical_model.fun_cons],
                            bounds=bounds_6D,
                            constraints_bounds=constraints_bounds_stdize_6D,
                            nb_it_bo=nb_it_cbo, 
                            normalization_data=normalization_data_6D)
        data_for_ns[-1].append((best_candidates, X_best_candidates))
            
        


# In[11]:


for itns,ns in enumerate(nb_samples):
    data_f_best_tmp = torch.vstack([x[0] for x in data_for_ns[itns]])
    data_f_best = list()
    for itbo in range(data_f_best_tmp.shape[1]):
        max_f = torch.max(data_f_best_tmp[:,itbo]).item()
        min_f = torch.min(data_f_best_tmp[:,itbo]).item()
        mean_f = torch.mean(data_f_best_tmp[:,itbo]).item()
        data_f_best.append(torch.tensor([min_f, mean_f, max_f], dtype=torch.float64))
    data_f_best = torch.hstack((torch.arange(1,data_f_best_tmp.shape[1]+1, dtype=torch.int).unsqueeze(-1),torch.vstack(data_f_best)))

    import pandas as pd
    df = pd.DataFrame(data_f_best.numpy()) #convert to a dataframe
    df.to_csv(f"CBO_6D_ns{ns:02d}.csv".format(ns),index=False) #save to file

