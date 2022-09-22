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

def calc_nnls_nut_aperp(ref_field_dir,ref_case_prefix_and_name):
    #tau = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_tau.npy'))
    S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_S.npy'))
    k = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_k.npy'))
    a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'))
    
    nut_nnls = np.empty(k.shape)
    r2_nnls = np.empty(k.shape)

    print('    Calculating nnls nut, aperp for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
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
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp11_nnls.npy'),aperp_nnls[:,0,0])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp12_nnls.npy'),aperp_nnls[:,0,1])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp13_nnls.npy'),aperp_nnls[:,0,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp22_nnls.npy'),aperp_nnls[:,1,1])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp23_nnls.npy'),aperp_nnls[:,1,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp33_nnls.npy'),aperp_nnls[:,2,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_r2_nnls.npy'),r2_nnls)

def calc_ls_nut_aperp(ref_field_dir,ref_case_prefix_and_name):
    #tau = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_tau.npy'))
    S = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_S.npy'))
    k = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_k.npy'))
    a = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'))
    
    nut_nnls = np.empty(k.shape)
    r2_nnls = np.empty(k.shape)

    print('    Calculating nnls nut, aperp for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    for i in range(len(S)):
        X = -2*np.delete(S[i,:].reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)
        y = np.delete(a[i,:].reshape((1,9)),[3,6,7],axis=1).flatten().reshape(-1,1)     
        reg_nnls = LinearRegression(positive=False,fit_intercept = False)
        reg_nnls.fit(X, y)
        nut_nnls[i] = reg_nnls.coef_[0,0]
        y_pred_nnls = reg_nnls.predict(X)
        r2_nnls[i] = r2_score(y, y_pred_nnls)
        if i % 10000 == 0:
            print('        Point: '+str(i)+'/'+str(len(S)))

    
    aperp_nnls = a + 2*np.repeat(nut_nnls,9).reshape(len(nut_nnls),3,3)*S 
    
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_nut_ls.npy'),nut_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp_ls.npy'),aperp_nnls)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp11_ls.npy'),aperp_nnls[:,0,0])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp12_ls.npy'),aperp_nnls[:,0,1])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp13_ls.npy'),aperp_nnls[:,0,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp22_ls.npy'),aperp_nnls[:,1,1])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp23_ls.npy'),aperp_nnls[:,1,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_aperp33_ls.npy'),aperp_nnls[:,2,2])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_r2_ls.npy'),r2_nnls)

    