#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

import os
import numpy as np
from scipy.interpolate import griddata
#from scipy.interpolate import RBFInterpolator
#from rbf.interpolate import RBFInterpolant

def interpolate_fields_fine_coarse(fine_field_dir,
                                          coarse_field_dir,
                                          output_dir,
                                          fine_case_prefix_and_name,
                                          coarse_case_prefix_and_name,
                                          output_case_prefix_and_name,
                                          fields_list,
                                          interp_method):
    print('[dataFoam] Interpolating fields from '+fine_field_dir+' to '+coarse_field_dir)
    C_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_C.npy'))

    C_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_C.npy'))
    
    print('[dataFoam] Fine field: '+fine_case_prefix_and_name+', '+str(len(C_coarse))+ ' points')
    print('[dataFoam] Coarse field: '+coarse_case_prefix_and_name+', '+str(len(C_coarse))+ ' points')

    for field_name in fields_list:
        print('[dataFoam] Interpolating '+field_name + ' using method '+interp_method)
        field_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_'+field_name+'.npy'))
        interp_field = griddata(C_fine,
                                field_fine,
                                C_coarse,
                                method=interp_method)

        if interp_method != 'nearest':
            ind_nan = np.argwhere(np.isnan(interp_field))
            if len(ind_nan) >0:
                print('[dataFoam] Interpolation found '+ str(len(ind_nan))+' nans, using nearest to fill these nans in')
            interp_field[ind_nan] = griddata(C_fine,
                                             field_fine,
                                             C_coarse[ind_nan], method='nearest')
        print('[dataFoam] Saving interpolated field to ' + str(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy')))
        np.save(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy'),interp_field)        
