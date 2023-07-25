#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
import multiprocessing as mp

def calc_k_b_a(ref_field_dir,ref_case_prefix_and_name):
    """Using the tau.npy field, calculates the TKE, anisotropy, and non-dimensional anisotropy tensor."""
    print('[dataFoam] Calculating k,b, and a for '+ref_case_prefix_and_name+ ' in '+ref_field_dir)
    tau = np.load(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_tau.npy'))
    k = 0.5*np.trace(tau,axis1=1,axis2=2)
    a = tau - 2/3 * k[:,None,None] * np.identity(3)
    b = a/(2*k[:,None,None])
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_k.npy'),k)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_a.npy'),a)
    np.save(os.path.join(ref_field_dir,ref_case_prefix_and_name+'_b.npy'),b)
    print('[dataFoam] Saved calculated results for k, a, b.')
    return
