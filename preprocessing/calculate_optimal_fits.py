#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
import multiprocessing as mp

def nnls_fit(S,a):
    """Calculates the non-negative least-squares fit for a = -2*nut*S. 
    Returns optimal nut, and the r2 value.
    """
    X = -2*np.delete(S.reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)
    y = np.delete(a.reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)  
    reg_nnls = LinearRegression(positive=True,fit_intercept = False)
    reg_nnls.fit(X, y)
    nut_nnls = reg_nnls.coef_[0,0]
    y_pred_nnls = reg_nnls.predict(X)
    r2_nnls = r2_score(y, y_pred_nnls)
    return nut_nnls, r2_nnls

def calc_nnls_nut_aperp(ref_field_dir,ref_case_prefix_and_name,parallel=True,LES_case=True):
    """Calls nnls_fit in parallel if needed, and also calculates aperp.
    ref: https://doi.org/10.1063/5.0083074
    """
    if LES_case:
        S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_SMean.npy'))
        a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aMean.npy'))
    else:
        S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_S.npy'))
        a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'))

    nut_nnls = np.empty(len(S))
    r2_nnls = np.empty(len(S))

    print('[dataFoam] Calculating nnls nut, aperp for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    if parallel:
        ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
        print(f'Parallel ncpus {ncpus}')
        pool = mp.Pool(processes=ncpus)
        results = pool.starmap(nnls_fit, [(Si,ai) for Si, ai in zip(S,a)])
        nut_nnls = np.array([result[0] for result in results])
        r2_nnls = np.array([result[1] for result in results])

    else:
        for i in range(len(S)):
            nut_nnls, r2_nnls = nnls_fit(S[i],a[i])
            if i % 10000 == 0:
                print('[dataFoam] Point: '+str(i)+'/'+str(len(S)))

    aperp_nnls = a + 2*np.repeat(nut_nnls,9).reshape(len(nut_nnls),3,3)*S  #a = -2*nut*S + aperp
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_nut_nnls.npy'),nut_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp_nnls.npy'),aperp_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_r2_nnls.npy'),r2_nnls)
    
    return 

def calc_modified_nnls_nut_aperp(ref_field_dir,ref_case_prefix_and_name,coarse_field_dir,coarse_case_prefix_and_name,parallel=True,LES_case=True):
    """Calculates aperp, using nut_rans as the optimal eddy viscosity."""
    if LES_case:
        S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_SMean.npy'))
        a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aMean.npy'))
    else:
        S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_S.npy'))
        a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'))
        k = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_k.npy'))
        k = np.maximum(k,np.ones(k.shape)*1E-10)

    nut_rans = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_nut.npy'))
    print('[dataFoam] Calculating MODIFIED nnls nut, aperp for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    aperp_nnls = a + 2*np.repeat(nut_rans,9).reshape(len(nut_rans),3,3)*S  #a = -2*nut*S + aperp
    bperp_nnls = aperp_nnls/(2*k[:,None,None])

    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_nut_nnls.npy'),nut_rans)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp_nnls.npy'),aperp_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_bperp_nnls.npy'),bperp_nnls)

    return
