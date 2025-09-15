import torch
import gpytorch
import botorch
import numpy as np
import pyDOE as doe
from loguru import logger

def default_values():
    out = {
        'nb_samples':5,
        'bounds': [[-0.2, -0.2], [0.2, 0.2]],
        'noise_lvl': 1e-4,
        'nb_it_optim': 100,
        'learning_rate': 1e-1,
        'nb_it_bo':20
    }
    return out

# various tools
class stdize:
    def __init__(self, input_data=None):
        self.internal_data = None
        if input_data is not None:
            self.output_data = self.init_stdize(input_data)          
    
    def init_stdize(self, x):
        input_shape = x.shape
        if len(input_shape) == 2:
            self.np = x.shape[1]
            x_tmp = x
        else:
            self.np = 1
            nb_item = input_shape[0]
            x_tmp = x.reshape(nb_item, 1)
        output_data = torch.zeros_like(x_tmp, dtype=torch.float64)
        self.internal_data = torch.zeros((2, self.np), dtype=torch.float64)
        for j in range(self.np):
            std, mean = torch.std_mean(x[:, j])
            output_data[:,j] = (x_tmp[:,j] - mean) / std
            self.internal_data[:,j] = torch.tensor([std, mean],dtype=torch.float64)
        if self.np == 1:
            output_data = output_data.reshape(input_shape)
        return output_data
    
    def stdize(self, x):
        if self.internal_data is None:
            raise ValueError("Standardization parameters not initialized.")
        input_shape = x.shape
        if self.np == 1:            
            x_tmp = x.reshape(np.prod(input_shape), 1)
        else:
            x_tmp = x
        output_data = torch.zeros_like(x_tmp, dtype=torch.float64)
        for j in range(self.np):
            mean = self.internal_data[1, j]
            std = self.internal_data[0, j]
            output_data[:,j] = (x_tmp[:,j]-mean)/std
        if self.np ==1:
            output_data = output_data.reshape(input_shape)
        return output_data

    def unstdize(self, x):
        if self.internal_data is None:
            raise ValueError("Standardization parameters not initialized.")
        input_shape = x.shape
        if self.np == 1:            
            x_tmp = x.reshape(np.prod(input_shape), 1)
        else:
            x_tmp = x
        output_data = torch.zeros_like(x_tmp, dtype=torch.float64)
        for j in range(self.np):
            mean = self.internal_data[1, j]
            std = self.internal_data[0, j]
            output_data[:,j] = std*x_tmp[:,j]+mean
        if self.np ==1:
            output_data = output_data.reshape(input_shape)
        return output_data
    
    def std_cons(self,constraints):
        pass
    
def lhs_distrib(bounds, nbs = default_values().get('nb_samples')):
    ndim = len(bounds[0])
    raw_samples = doe.lhs(ndim, samples=nbs, criterion='maximin', iterations=1000)
    # scale samples
    sample_points = torch.zeros((nbs, ndim), dtype=torch.float64)
    for i in range(ndim):
        sample_points[:, i] = torch.tensor(raw_samples[:, i] * (bounds[1][i] - bounds[0][i]) + bounds[0][i], dtype=torch.float64)
    return sample_points

def init_GP_structure():
    likelihood_obj = gpytorch.likelihoods.FixedNoiseGaussianLikelihood
    mean_obj = gpytorch.means.LinearMean #gpytorch.means.ConstantMean
    covar_obj = gpytorch.kernels.MaternKernel
    GP_type_obj = botorch.models.SingleTaskGP
    mll_obj = gpytorch.mlls.ExactMarginalLogLikelihood
    optim_obj = torch.optim.Adam
    return likelihood_obj, mean_obj, covar_obj, GP_type_obj, mll_obj, optim_obj

# initialize surrogate model
def init_surrogate_model(X,Z, 
                         nb_it_optim=default_values().get('nb_it_optim'), 
                         noise_lvl=default_values().get('noise_lvl'), 
                         learning_rate=default_values().get('learning_rate'),
                         verbose=1):
    likelihood_obj, mean_obj, covar_obj, GP_type_obj, mll_obj, optim_obj = init_GP_structure()

    # batch size
    batch_shape = torch.Size([Z.shape[0]])
    # likelihood
    likelihood_module = likelihood_obj(
        noise = torch.ones(X.shape[0]) * noise_lvl,
        num_tasks = 1,
        rank = 0,
        #    batch_shape = batch_shape
        )
    # mean
    mean_module = mean_obj(input_size = X.shape[1],
                           #batch_shape = batch_shape,
                           bias = False)
    # covariance
    covariance_module = covar_obj()#batch_shape = batch_shape)
    # build GP model
    gp = GP_type_obj(X, Z, 
                     likelihood = likelihood_module,
                     covar_module = covariance_module,
                     mean_module = mean_module,
                     outcome_transform=None,
                     input_transform=None)
    # likelihood evaluation
    mll = mll_obj(gp.likelihood, gp)
    # botorch.fit.fit_gpytorch_mll(mll)
    mll_train = mll.to(X)
    
    # optimizer
    optimizer_module = optim_obj(gp.parameters(), 
                                 lr = learning_rate)
    
    # start training
    gp.train()
    for it in range(nb_it_optim):
        # clear gradients
        optimizer_module.zero_grad()
        # forward pass through the model to obtain output MulivariateNormal
        out = gp(gp.train_inputs[0])
        # compute negative marginal log likelihood
        loss = -mll_train(out, gp.train_targets).sum()
        # back propagation for gradients
        loss.backward()
        optimizer_module.step()
        if verbose >1:
            logger.info('Iter %d/%d - Loss: %.3f   lengthscale: %.3f   noise: %.3f' % (
        it + 1, nb_it_optim, loss.item(),
        covariance_module.lengthscale.item(),
        likelihood_module.noise.max()
    ))

    gp.eval()
        
    return gp


def find_best_candidate(X, Zval, Cval=None, constraints_bounds= None):
    if Cval is None:
        violated = [torch.tensor([False]*len(Zval))]
    else:
        violated = [torch.tensor([False]*len(Zval))]*len(Cval)
    if constraints_bounds:
        # check that constraints has been satisfied
        # check all Cval        
        for cons_it,cons_i_val in enumerate(Cval):
            for itC,iCval in enumerate(cons_i_val):
                current_bounds = constraints_bounds[cons_it]
                if (iCval < current_bounds[0]) or (iCval > current_bounds[1]):
                    violated[cons_it][itC]=True
    # concatenate violations
    violated_any = torch.vstack(violated).any(dim=0)
    #
    IXs = torch.argsort(Zval.flatten())
    violated_sorted = violated_any[IXs]
    for (v,iiX) in zip(violated_sorted, IXs):
        if not v:
            IX = iiX
            break
    best_candidate = Zval[IX]
    Xbest = X[IX, :]
    return best_candidate,IX, Xbest

def execute_CBO(gp_obj_current,                
                fun_obj,
                fun_cons=None,
                bounds=None,
                constraints_bounds=None, 
                nb_it_bo=default_values().get('nb_it_bo'),
                verbose=1,
                XY_plot=None,
                normalization_data=None):
    # constraints
    cbo_activated = fun_cons is not None and constraints_bounds is not None
    nb_constraints = 0
    if fun_cons is not None:
        if isinstance(fun_cons, str):
            nb_constraints = 1
        else:
            nb_constraints = len(fun_cons)
    #load normalization data
    dataStdizeX = normalization_data.get('X')
    dataStdizeZ = normalization_data.get('Z')
    dataStdizeC = normalization_data.get('C')
    #
    bounds_BO = torch.tensor(bounds, dtype=torch.float64)
    GP_obj_at_grid = list()
    acq_values_opt = list()
    best_candidates = list()
    X_best_candidates = list()
    if cbo_activated:
        GP_cons_at_grid = [list()]*len(fun_cons)
    #
    acq_at_grid = list()
    for it in range(nb_it_bo):
        logger.info(f"Iteration {it+1}/{nb_it_bo}")
        #get training data
        if nb_constraints > 0:
            X = gp_obj_current.train_inputs[0][0]
        else:
            X = gp_obj_current.train_inputs[0]
        ZC_multi_raw = gp_obj_current.train_targets
        if hasattr(gp_obj_current,'outcome_transform'):
            unraw, _ = gp_obj_current.outcome_transform.untransform(Zraw)
            ZC_multi = ZC_multi_raw[0]
        else:
            ZC_multi = ZC_multi_raw
        #
        C = None
        if cbo_activated:
            Z = ZC_multi[0]
            C = ZC_multi[1:(nb_constraints+1)]
        else:
            Z= ZC_multi
        # get best candidate
        best_candidate,_,X_best = find_best_candidate(X, Z, Cval=C, constraints_bounds=constraints_bounds)
        best_candidates.append(dataStdizeZ.unstdize(best_candidate.unsqueeze(0).detach()))
        X_best_candidates.append(dataStdizeX.unstdize(X_best.unsqueeze(0).detach()))
        
        # initialize acquisition function
        if constraints_bounds is None:
            acq_func = botorch.acquisition.ExpectedImprovement(gp_obj_current,
                                                                  best_f=best_candidate,
                                                                  maximize=False)
        else:
            # adaptation of constraints_bounds
            formatted_cons_bnds = dict()
            for ibnd,vbnd in enumerate(constraints_bounds):
                if vbnd[0].abs() == torch.inf:
                    lb = None
                else:
                    lb = vbnd[0]
                if vbnd[1].abs() == torch.inf:
                    ub = None
                else:
                    ub = vbnd[1]
                formatted_cons_bnds.update({ibnd+1: (lb,ub)})

            # botorch.fit_gpytorch_mll(mll)
            acq_func = botorch.acquisition.ConstrainedExpectedImprovement(gp_obj_current,
                                                                          objective_index=0,
                                                                          best_f=best_candidate,
                                                                          constraints=formatted_cons_bnds,
                                                                          maximize=False)  
        # run optimization
        bounds_acqf=dataStdizeX.stdize(bounds_BO)
        candidate, eimax = botorch.optim.optimize_acqf(acq_function=acq_func,      
                                                      q=1,
                                                      bounds=bounds_acqf,
                                                      num_restarts=20,
                                                      raw_samples=200)
                                                    #   ,options={})
        acq_values_opt.append(eimax)
        # evaluation of the actual function
        new_X = candidate.detach()
        new_X_un = dataStdizeX.unstdize(new_X)    
        # filter/enforce new sample point in bounds
        new_X_un = torch.clamp(new_X_un, min=bounds_BO[0], max=bounds_BO[1])
        new_Z_tmp = fun_obj(new_X_un)
        if len(new_Z_tmp.shape)==0:
            new_Z_tmp = new_Z_tmp.reshape(1)
        new_Z_un = torch.tensor(new_Z_tmp, dtype=torch.float64)
        new_Z = dataStdizeZ.stdize(new_Z_un)    
        if cbo_activated:
            new_C_un = list()
            for fun_cons_i in fun_cons:
                new_C_tmp = fun_cons_i(new_X_un)
                if len(new_C_tmp.shape)==0:
                    new_C_tmp = new_C_tmp.reshape(1)
                new_C_un.append(torch.tensor(new_C_tmp, dtype=torch.float64))
            new_C_un = torch.tensor(new_C_un, dtype=torch.float64)
            new_C = dataStdizeC.stdize(new_C_un) 
        if verbose>0:
            txt = f"New sample point {new_X_un.numpy()} with objective value {new_Z_un.numpy()}"
            if cbo_activated:
                txt += " and constraints values " + ", ".join([f"c{i}: {c.numpy()}" for i, c in enumerate(new_C_un)])
            logger.info(txt)
        
        # 
        X_update = torch.vstack([X, new_X])    
        Z_update = torch.cat([Z, new_Z]).unsqueeze(-1)
        ZC_update = Z_update
        if cbo_activated:
            C_update = torch.hstack((C,new_C.unsqueeze(-1))).unsqueeze(-1)
            ZC_update = torch.hstack((Z_update,*C_update))
        # update GP
        gp_obj_current = init_surrogate_model(X_update.double(), 
                                            ZC_update.double(),
                                            noise_lvl=1e-4, verbose=verbose)
        
        # evaluate functions and GPs at provided grid
        if XY_plot is not None:
            ZC_multi_GP = gp_obj_current(XY_plot).mean.detach()            
            if cbo_activated:
                GP_obj_at_grid.append(ZC_multi_GP[0])
                for ic in range(nb_constraints):
                    GP_cons_at_grid[ic].append(ZC_multi_GP[ic+1])
            else:
                GP_obj_at_grid.append(ZC_multi_GP)
            Xtmp = XY_plot.reshape(XY_plot.shape[0],1,XY_plot.shape[1])
            acq_at_grid.append(acq_func(Xtmp).detach())
        
        # check interpolation
        ZC_GP_at_samples = gp_obj_current(X_update).mean.detach()
        if cbo_activated:
            Z_GP_at_samples = ZC_GP_at_samples[0]
        else:
            Z_GP_at_samples = ZC_GP_at_samples
        logger.info(f'Z maxi interpolation error {(Z_update.flatten()-Z_GP_at_samples).double().abs().max()}')
        if cbo_activated:
            for ic in range(nb_constraints):
                C_GP_at_samples = ZC_GP_at_samples[ic+1]
                logger.info(f'{ic} C maxi interpolation error {(ZC_update[:,ic+1].flatten()-C_GP_at_samples).double().abs().max()}')
    if not cbo_activated:
        return X_update, \
            Z_update, \
                gp_obj_current, \
                    GP_obj_at_grid, \
                        acq_at_grid, \
                            torch.tensor(acq_values_opt, dtype=torch.float64), \
                                torch.tensor(best_candidates, dtype=torch.float64), \
                                    torch.vstack(X_best_candidates)
    else:
        return X_update,\
            ZC_update,\
                gp_obj_current,\
                    GP_obj_at_grid,\
                        GP_cons_at_grid,\
                            acq_at_grid,\
                                torch.tensor(acq_values_opt, dtype=torch.float64),\
                                    torch.tensor(best_candidates, dtype=torch.float64),\
                                        torch.vstack(X_best_candidates)