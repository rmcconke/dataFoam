#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

import os
import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
import multiprocessing as mp
import time
def calc_k_b_a(ref_field_dir,ref_case_prefix_and_name):
    print('    Calculating k,b, and a for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    tau = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_tau.npy'))
    k = 0.5*np.trace(tau,axis1=1,axis2=2)
    a = tau - 2/3 * k[:,None,None] * np.identity(3)
    b = a/(2*k[:,None,None])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_k.npy'),k)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'),a)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_b.npy'),b)
    print('    Saved calculated results for k, a, b.')

def calc_S_R(ref_field_dir,ref_case_prefix_and_name):
    print('    Calculating S and R for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    gradU = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_gradU.npy'))
    S = 0.5*(gradU + np.transpose(gradU,(0,2,1)))
    R = -0.5*(gradU - np.transpose(gradU,(0,2,1))) #negative sign because OpenFOAM outputs the transpose of the Jacobian. negative sign matches convention used in RANS callculations.
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_S.npy'),S)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_R.npy'),R)
    print('    Saved calculated results for S and R.')

def nnls_fit(S,a):
    X = -2*np.delete(S.reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)
    y = np.delete(a.reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)  
    reg_nnls = LinearRegression(positive=True,fit_intercept = False)
    reg_nnls.fit(X, y)
    nut_nnls = reg_nnls.coef_[0,0]
    y_pred_nnls = reg_nnls.predict(X)
    r2_nnls = r2_score(y, y_pred_nnls)
    return nut_nnls, r2_nnls

def calc_nnls_nut_aperp(ref_field_dir,ref_case_prefix_and_name,parallel=True):
    #tau = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_tau.npy'))
    S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_SMean.npy'))
    a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aMean.npy'))

    nut_nnls = np.empty(len(S))
    r2_nnls = np.empty(len(S))

    print('    Calculating nnls nut, aperp for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    if parallel:
        ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
        print(f'Parallel ncpus {ncpus}')
        pool = mp.Pool(processes=ncpus)
        results = pool.starmap(nnls_fit, [(Si,ai) for Si, ai in zip(S,a)])
        nut_nnls = np.array([result[0] for result in results])
        r2_nnls = np.array([result[1] for result in results])

    else:
        for i in range(len(S)):
            X = -2*np.delete(S[i,:].reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)
            y = np.delete(a[i,:].reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)     
            reg_nnls = LinearRegression(positive=True,fit_intercept = False)
            reg_nnls.fit(X, y)
            nut_nnls[i] = reg_nnls.coef_[0,0]
            y_pred_nnls = reg_nnls.predict(X)
            r2_nnls[i] = r2_score(y, y_pred_nnls)
            if i % 10000 == 0:
                print('        Point: '+str(i)+'/'+str(len(S)))

    aperp_nnls = a + 2*np.repeat(nut_nnls,9).reshape(len(nut_nnls),3,3)*S 
    
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_nut_nnls.npy'),nut_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp_nnls.npy'),aperp_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_r2_nnls.npy'),r2_nnls)
    mag_aperp_nnls = np.linalg.norm(aperp_nnls,axis=(1,2))
    print(f'Magnitude of a_perp shape: {mag_aperp_nnls.shape}')
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp_nnls_magnitude.npy'),mag_aperp_nnls)

def gn_fit(yi, T_colsi, N, alpha):
    Xi = np.empty((6,N-1))
    for row in range(6):
        Xi[row,:] = T_colsi[:,row]
    reg = Ridge(alpha=alpha,fit_intercept = False)
    reg.fit(Xi, yi)
    #print(f'Coefficient shape: {reg.coef_.shape}')
    gn = reg.coef_
    gn_error = yi - reg.predict(Xi)
    gn_b_pred = reg.predict(Xi)
    return gn, gn_error, gn_b_pred

def calc_gn(N, alpha, mapped_field_dir, mapped_case_prefix_and_name, coarse_field_dir, coarse_case_prefix_and_name,parallel=True):
    # Load aperp
    aperp = np.load(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+'_aperp_nnls.npy'))
    k = np.load(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+'_kMean_tauMean.npy'))
    # Flatten/reshape aperp
    bperp = aperp/(2*k[:,None,None])
    y = np.column_stack((bperp[:,0,0],bperp[:,0,1],bperp[:,0,2],bperp[:,1,1],bperp[:,1,2],bperp[:,2,2]))

    # Load T's
    # Flatten/reshape T's
    print('    Calculating gn fit for '+mapped_case_prefix_and_name+ ' and '+coarse_case_prefix_and_name)

    T_cols = np.empty((len(y),N-1,6))
    for n in range(0,N-1):
        Tn = np.load(os.path.join(coarse_field_dir,f'{coarse_case_prefix_and_name}_T{n+2}.npy'))
        Xn = np.column_stack((Tn[:,0,0],Tn[:,0,1],Tn[:,0,2],Tn[:,1,1],Tn[:,1,2],Tn[:,2,2]))
        T_cols[:,n,:]=Xn

    gn = np.empty((len(y),N-1))
    gn_error = np.empty((len(y),6))
    gn_b_pred = np.empty((len(y),6))

    if parallel:
        ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
        print(f'Parallel ncpus {ncpus}')
        pool = mp.Pool(processes=ncpus)
        results = pool.starmap(gn_fit, [(yi,T_colsi, N, alpha) for yi, T_colsi in zip(y,T_cols)])
        gn = np.array([result[0] for result in results])
        gn_error = np.array([result[1] for result in results])
        gn_b_pred = np.array([result[2] for result in results])

    for i in range(10):#range(len(y)):
        yi = y[i,:]
        Xi = np.empty((6,N-1))
        for row in range(6):
            Xi[row,:] = T_cols[i,:,row]
        reg = Ridge(alpha=alpha,fit_intercept = False)
        reg.fit(Xi, yi)
        #print(f'Coefficient shape: {reg.coef_.shape}')
        # Calculate errors
        gn[i,:] = reg.coef_
        gn_b_pred[i,:] = reg.predict(Xi)
        gn_error[i,:] = yi - reg.predict(Xi)
        #print(f'yi: {yi}')
        #print(f'pd: {reg.predict(Xi)}')
        #print(f'gn: {gn[i,:]}')
        #print(f'ge: {gn_error[i,:]}')

    gn_error_tensor = np.stack((gn_error[:,0], gn_error[:,1], gn_error[:,2], gn_error[:,1], gn_error[:,3], gn_error[:,4], gn_error[:,2], gn_error[:,4],gn_error[:,5]),axis=1)
    gn_error_tensor = gn_error_tensor.reshape((len(gn_error_tensor),3,3))

    gn_b_pred_tensor = np.stack((gn_b_pred[:,0], gn_b_pred[:,1], gn_b_pred[:,2], gn_b_pred[:,1], gn_b_pred[:,3], gn_b_pred[:,4], gn_b_pred[:,2], gn_b_pred[:,4],gn_b_pred[:,5]),axis=1)
    gn_b_pred_tensor = gn_b_pred_tensor.reshape((len(gn_error_tensor),3,3))
    gn_a_pred_tensor = gn_b_pred_tensor *(2*k[:,None,None])

    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}.npy'),gn)
    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}_error.npy'),gn_error_tensor)
    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}_error_magnitude.npy'),np.linalg.norm(gn_error_tensor,axis=(1,2)))
    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}_error_relative_magnitude.npy'),np.linalg.norm(gn_error_tensor,axis=(1,2))/np.maximum(1E-8,np.linalg.norm(bperp,axis=(1,2))))

    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}_b_pred.npy'),gn_b_pred_tensor)
    np.save(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gn_N{N}_alpha{alpha}_a_pred.npy'),gn_a_pred_tensor)


def calc_gamma(mapped_field_dir, mapped_case_prefix_and_name, coarse_field_dir, coarse_case_prefix_and_name):
    nut_rans=np.load(os.path.join(coarse_field_dir,f'{coarse_case_prefix_and_name}_nut.npy'))
    nut_les=np.load(os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_nut_nnls.npy'))
    nut_rans_divide = np.maximum(1E-8,nut_rans)
    print(np.amin(nut_rans_divide))
    nut_les_divide = np.maximum(1E-8,nut_les)
    print(np.amin(nut_les_divide))
    gamma=np.log(nut_les_divide/nut_rans_divide)
    gamma_save_path = os.path.join(mapped_field_dir,mapped_case_prefix_and_name+f'_gamma.npy')
    np.save(gamma_save_path,gamma)
    print(f'Saved gamma to {gamma_save_path}.')
