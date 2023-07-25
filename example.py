
from utilities.MLDatasetFromFoamCase import MLDatasetFromFoamCase
import os
import numpy as np
from preprocessing.assemble_dataframe import assemble_dataframe

# Simple example of converting an OpenFOAM case into a collection of numpy binaries, and then a csv file.
dataFoam_folder = os.getcwd()

# This code writes numpy binaries from an OpenFOAM case. 
# The fields are written in OpenFOAM using the custom utilies in foam_applications.
# The OpenFOAM utlility called depends on the write_fields_application - different utilities are needed depending on if the data comes from RANS, LES, or DNS.
# The list of fields read and saved for each case_type are given in the foam_field_list for a given case_type in MLDatasetFromFoamCase.
foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'example_data/numpy/komegasst'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'example_data'),
                                       case_name='case_1p0',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0') # This saves the numpy binaries to MLDatasetFromFoamCase.data_save_path

foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'example_data/numpy/DNS'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'example_data'),
                                       case_name='case_1p0_DNS',
                                       case_type='DNS',
                                       write_fields_application='writeFields_DNS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='DNS_case_1p0') # This saves the numpy binaries to MLDatasetFromFoamCase.data_save_path
# If your fine data is not available as an OpenFOAM case, you can write these biniaries yourself using another code, and then continue

# After writing a collection of numpy binaries, we often need to interpolate the fine field onto the coarse field.
# The interpolate_fields_fine_coarse function interpolates the fine fields in fields_list onto the coordinates given by "{coarse_case_prefix_and_name}_C.npy"
# The coarse and fine coordinates are automatically written in MLDatasetFromFoamCase 
from preprocessing.map_coarse_fine import interpolate_fields_fine_coarse

# This code saves the interpolated fields with name {output_case_prefix_and_name}_{field_name}.npy to output_dir
interpolate_fields_fine_coarse(fine_field_dir = os.path.join(dataFoam_folder,'example_data/numpy/DNS'),
                               coarse_field_dir = os.path.join(dataFoam_folder,'example_data/numpy/komegasst'),
                               output_dir = os.path.join(dataFoam_folder,'example_data/numpy/DNS_mapped'),
                               fine_case_prefix_and_name='DNS_case_1p0',
                               coarse_case_prefix_and_name='komegasst_case_1p0',
                               output_case_prefix_and_name='DNS_mapped_case_1p0',
                               fields_list=['U','k','TauDNS'], 
                               interp_method = 'nearest' # Usually, linear is recommended. For the example, the DNS data and RANS data are already on the same grid.
                              )

assemble_dataframe(data_parent_dir=os.path.join(dataFoam_folder,'example_data/numpy'),
                   field_dict = {'komegasst': ['U','C','nut','k','S'],
                                 'DNS_mapped': ['U','k','TauDNS']},
                  case_name_list = ['case_1p0'], # Additional cases can be added to this csv file, as most studies use multiple cases
                  output_file = os.path.join(dataFoam_folder,'example_data/output_example.csv'))








