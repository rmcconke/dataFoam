#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
def assemble_dataframe(data_parent_dir: str,
                       field_dict: dict,
                       case_name_list: list, 
                       output_file: str):
    """
    Generates a csv file from a pandas dataframe, which in turn is an assembly of various numpy fields.
    data_parent_dir: parent directory where numpy fields are stored. Subdirectories in this directory should have the generation method as their name, e.g.
    - data_parent_dir
        - komegasst
        - DNS
        - LES
    field_dict: a dictionary to determine which fields should be stored in the csv file for each generation method.
    field_dict example:
    field_dict = {'komegasst': ['omega'],
                      'nnls_interp':['a_perp_nnls','nut_nnls']}
    case_name_list: list of flow cases to include in the csv file.
    output_file: output filename for the csv file. Should end in .csv.
    """

    print('[dataFoam] Assembling dataframe for cases: '+str(case_name_list)+', outputting final csv file to '+output_file)

    dataframe_main = pd.DataFrame()

    # Loop over cases
    for case in case_name_list:
        
        print(f'    Case: {case}')

        # dataframe containing data from current case
        dataframe_case = pd.DataFrame()

        # Loop over field generation method
        for field_dir in field_dict:
            #
            print(f'        {field_dir}')

            # Loop over field name
            for field_name in field_dict[field_dir]:
                print(f'            {field_name}')

                # Try to load the numpy field, raise fatal error if it can't be found
                try: 
                    field = np.float32(np.load(os.path.join(data_parent_dir,field_dir,field_dir+'_'+case+'_'+field_name+'.npy')))
                except:
                    raise LookupError('Cannot find field '+field_name)
                
                field_type = get_field_type(field_name,field)

                # Depending on the field type, add different entry format to the dataframe. Since the dataframe is tabular, it needs columns, and tensors need to be split up into individual entries.
                
                # Scalars can be added directly
                if field_type == 'scalar':
                    dataframe_case = add_data(dataframe_case,field_dir,field_name,field)

                # Vectors and matrices should be added one component at a time
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

                # The Tensor basis field requires special handling, since it is of shape [N,10,3,3]
                if field_type == 'tensor_basis':
                    for i in range(10):
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_11',field[:,i,0,0])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_12',field[:,i,0,1])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_13',field[:,i,0,2])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_22',field[:,i,1,1])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_23',field[:,i,1,2])
                        dataframe_case = add_data(dataframe_case,field_dir,field_name+'_'+str(i+1)+'_33',field[:,i,2,2])

        # Convert to float 32
        dataframe_case=dataframe_case.astype('float32')

        # Add a Case entry to the dataframe
        dataframe_case['Case'] = case

        # Concatenate dataframe into the main dataframe
        dataframe_main = pd.concat((dataframe_main,dataframe_case),axis=0)

    dataframe_main['Case']=dataframe_main['Case'].astype('string')
    print('Final dataframe columns:')
    print(dataframe_main.columns)
    print('Memory info for final dataframe: ')
    print(dataframe_main.memory_usage(deep=True).to_string())
    print(str(dataframe_main.dtypes.to_string()))
    print('Saving dataframe to '+output_file)
    dataframe_main.to_csv(output_file)
    print('Finished saving dataframe.')

def get_field_type(fieldname,field):
    # Determine field type based on number of dimensions, and for symmetric matrices, the fieldname itself.
    ndims = field.ndim
    if ndims == 1:
        fieldtype = 'scalar'
    if ndims == 2: 
        fieldtype = 'vector'
    if ndims == 3:
        if fieldname in ['S','Shat','b','tau','a']:
            fieldtype = 'symm_matrix'
        else: 
            fieldtype = 'matrix'
    if ndims == 4 and field.shape[1]==10:
        fieldtype = 'tensor_basis'

    return fieldtype

def add_data(dataframe_case,field_dir,field_name,field):
    # Expects field to be a 1D scalar
    # Usually throws performance warnings after several calls, but performance is usually acceptable.
    dataframe_case[field_dir+'_'+field_name] = field.tolist()
    return dataframe_case