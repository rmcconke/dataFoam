#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from scipy.interpolate import griddata

def interpolate_fields_fine_coarse(fine_field_dir,
                                          coarse_field_dir,
                                          output_dir,
                                          fine_case_prefix_and_name,
                                          coarse_case_prefix_and_name,
                                          output_case_prefix_and_name,
                                          fields_list,
                                          interp_method):
    """Interpolate an LES/DNS field onto a RANS field.
    *_field_dir: contains numpy files e.g. ..../DNS/ ofr ..../komegasst/
    output_dir: directory for saving interpolated fields
    *_case_prefix_and_name: prefix (e.g. komegasst) and case name. e.g. komegasst_case_1p0 or DNS_case_1p0
    fields_list: list of DNS/LES fields to interpolate 
    interp_method: method used by scipy.interpolate.griddata. Suggested to use linear or nearest. The code will fill any NaNs with nearest.
    NaNs may be caused by asking certain interpolation methods to extrapolate.
    """

    print('[dataFoam] Interpolating fields from '+fine_field_dir+' to '+coarse_field_dir)
    C_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_C.npy'))
    C_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_C.npy'))

    print('[dataFoam] Fine field: '+fine_case_prefix_and_name+', '+str(len(C_fine))+ ' points')
    print('[dataFoam] Coarse field: '+coarse_case_prefix_and_name+', '+str(len(C_coarse))+ ' points')

    for field_name in fields_list:
        print('[dataFoam] Interpolating '+field_name + ' using method '+interp_method)
        field_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_'+field_name+'.npy'))
        interp_field = griddata(C_fine,
                                field_fine,
                                C_coarse,
                                method=interp_method)
        # Below code fills NaNs with nearest, and prints info about how many NaNs were found
        if interp_method != 'nearest':
            ind_nan = np.argwhere(np.isnan(interp_field))
            if len(ind_nan) >0:
                print('[dataFoam] Interpolation found '+ str(len(ind_nan))+' nans, using nearest to fill these nans in')
            interp_field[ind_nan] = griddata(C_fine,
                                             field_fine,
                                             C_coarse[ind_nan], method='nearest')
        print('[dataFoam] Saving interpolated field to ' + str(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy')))
        np.save(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy'),interp_field)    
    return    

def interpolate_fields_fine_coarse_DUCT(fine_field_dir,
                                          coarse_field_dir,
                                          output_dir,
                                          fine_case_prefix_and_name,
                                          coarse_case_prefix_and_name,
                                          output_case_prefix_and_name,
                                          fields_list,
                                          interp_method):
    """
    Special version of interpolate_fields_fine_coarse for a duct case.
    Given the current order of field processing, the duct cases require interpolation on the 2D yz plane, which is achieved by setting C_fine_x to zero, and doing 2D interpolation in griddata. 
    This problem arises from the two data-containing yz planes not having the same x coordinate.
    """
    print('[dataFoam] Interpolating DUCT fields from '+fine_field_dir+' to '+coarse_field_dir)
    C_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_C.npy'))
    C_fine[:,0] = np.zeros(len(C_fine))
    C_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_C.npy'))
    
    print('[dataFoam] Fine field: '+fine_case_prefix_and_name+', '+str(len(C_fine))+ ' points')
    print('[dataFoam] Coarse field: '+coarse_case_prefix_and_name+', '+str(len(C_coarse))+ ' points')

    for field_name in fields_list:
        print('[dataFoam] Interpolating '+field_name + ' using method '+interp_method)
        field_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_'+field_name+'.npy'))
        interp_field = griddata(C_fine[:,1:],
                                field_fine,
                                C_coarse[:,1:],
                                method=interp_method)

        if interp_method != 'nearest':
            ind_nan = np.argwhere(np.isnan(interp_field))
            if len(ind_nan) >0:
                print('[dataFoam] Interpolation found '+ str(len(ind_nan))+' nans, using nearest to fill these nans in')
                interp_field[ind_nan] = griddata(C_fine[:,1:],
                                                field_fine,
                                                C_coarse[ind_nan][:,1:], method='nearest')
        print('[dataFoam] Saving interpolated field to ' + str(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy')))
        np.save(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy'),interp_field)        
