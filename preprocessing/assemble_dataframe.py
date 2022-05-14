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
import pandas as pd
def assemble_dataframe(data_parent_dir,
                       field_dict,
                       case_name_list, 
                       output_file):
    #field_dict example (only scalars work at this time):
    #field_dict = {'komegasst': ['omega'],
    #                  'nnls_interp':['a_perp_nnls','nut_nnls']}
    print('Assembling dataframe for cases: '+str(case_name_list)+', outputting to '+output_file)

    dataframe_main = pd.DataFrame()

    for case in case_name_list:
        print(case)
        dataframe_case = pd.DataFrame()
        for field_dir in field_dict:
            print('    '+field_dir)
            for field_name in field_dict[field_dir]:
                print('        '+field_name)
                try: 
                    field = np.load(os.path.join(data_parent_dir,field_dir,field_dir+'_'+case+'_'+field_name+'.npy'))
                except:
                    raise LookupError('Cannot find field '+field_name)
                field_type = get_field_type(field_name,field)
                if field_type == 'scalar':
                    dataframe_case = add_data(dataframe_case,field_dir,field_name,field)
                if field_type == 'vector':
                    for i in range(field.shape[1]):
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1),field[:,i])
                if field_type == 'matrix':
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_11',field[:,0,0])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_12',field[:,0,1])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_13',field[:,0,2])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_21',field[:,1,0])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_22',field[:,1,1])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_23',field[:,1,2])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_31',field[:,2,0])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_32',field[:,2,1])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_33',field[:,2,2])
                    
                if field_type == 'symm_matrix':
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_11',field[:,0,0])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_12',field[:,0,1])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_13',field[:,0,2])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_22',field[:,1,1])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_23',field[:,1,2])
                    dataframe_case = add_data(dataframe_case,field_dir,field_name+'_33',field[:,2,2])

                if field_type == 'tensor_basis':
                    for i in range(10):
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_11',field[:,i,0,0])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_12',field[:,i,0,1])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_13',field[:,i,0,2])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_22',field[:,i,1,1])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_23',field[:,i,1,2])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_33',field[:,i,2,2])
            
        dataframe_case['Case'] = case
        dataframe_main = pd.concat((dataframe_main,dataframe_case),axis=0)
    print('    Saving dataframe to '+output_file)
    dataframe_main.to_csv(output_file)

def get_field_type(fieldname,field):
    ndims = field.ndim
    if ndims == 1:
        fieldtype = 'scalar'
    if ndims == 2: 
        fieldtype = 'vector'
    if ndims == 3:
        if fieldname in ['S','Shat','b','tau']:
            fieldtype = 'symm_matrix'
        else: 
            fieldtype = 'matrix'
    if ndims == 4 and field.shape[1]==10:
        fieldtype = 'tensor_basis'

    return fieldtype

def add_data(dataframe_case,field_dir,field_name,field):
    # Expects field to be a 1D scalar
    dataframe_case[field_dir+'_'+field_name] = field.tolist()
    # print('    Added '+field_dir+'_'+field_name+' into dataframe, shape: '+str(field.shape))
    return dataframe_case