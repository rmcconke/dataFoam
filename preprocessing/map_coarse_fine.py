#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

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
    print('    Interpolating fields from '+fine_field_dir+' to '+coarse_field_dir)
    x_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_Cx.npy'))
    y_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_Cy.npy'))
    z_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_Cz.npy'))
    
    x_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_Cx.npy'))
    y_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_Cy.npy'))
    z_coarse = np.load(os.path.join(coarse_field_dir,coarse_case_prefix_and_name+'_Cz.npy'))
    
    print('    Fine field: '+fine_case_prefix_and_name+', '+str(len(x_fine))+ ' points')
    print('    Coarse field: '+coarse_case_prefix_and_name+', '+str(len(x_coarse))+ ' points')

    for field_name in fields_list:
        print('    Interpolating '+field_name + ' using method '+interp_method)
        field_fine = np.load(os.path.join(fine_field_dir,fine_case_prefix_and_name+'_'+field_name+'.npy'))
        interp_field = griddata((x_fine,y_fine,z_fine),
                                field_fine,
                                (x_coarse,y_coarse,z_coarse),
                                method=interp_method)
        
        if interp_method != 'nearest':
            ind_nan = np.argwhere(np.isnan(interp_field))
            if len(ind_nan) >0:
                print('    Interpolation found '+ str(len(ind_nan))+' nans, using nearest to fill these nans in')
            interp_field[ind_nan] = griddata((x_fine,y_fine,z_fine),
                                             field_fine,
                                             (x_coarse[ind_nan],y_coarse[ind_nan]), method='nearest')
        print('    Saving interpolated field to ' + str(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy')))
        np.save(os.path.join(output_dir,output_case_prefix_and_name+'_'+field_name+'.npy'),interp_field)        
