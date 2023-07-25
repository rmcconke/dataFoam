#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

import os
import numpy as np
from utilities.foamIO.readFoam import readFoamField, get_endtime, get_cell_centres
from utilities.foamIO.changeFoamSystemDictEntry import changeFoamSystemDictEntry

class MLDatasetFromFoamCase: 
    """ 
    Class for a foam case that we want to store a numpy dataset for. 
    The fields stored depend on the case_type. the self.foam_field_list attribute is used to specify the fields which are converted from foam to numpy.
    RANS case types store the full set of invariants, while reference cases (LES/DNS) only store U/gradU/tau related fields.
    """

    def __init__(self,data_save_path,foam_parent_dir,case_name,case_type,write_fields_application,write_fields_flag=True):
        """Constructor
        data_save_path: output numpy folder
        foam_parent_directory - e.g., the komegasst foam dataset directory, which contains all of the komegasst cases
        case_name - unique identifier, e.g. case_1p0
        case_type: kepsilonphitf, komegasst, reference (see if statements below)
        write_fields_application: name of OpenFOAM application called to write the additional fields (see dataFoam/foam_applications)
        overwrite_flag: whether to call the write_fields_application, which may be time consuming
        """

        print('[dataFoam] Initializing MLDatasetFromFoamCase....')
        self.case_type = case_type
        self.foam_parent_dir = foam_parent_dir
        self.directory = os.path.join(self.foam_parent_dir,case_name)
        print('[dataFoam] Case directory: '+self.directory)
        self.writeFieldsApplication=write_fields_application
        self.case_name = case_name
        self.data_save_path = data_save_path
        self.write_fields_flag=write_fields_flag

        self.endtime = get_endtime(self.directory)

        
        if self.case_type == 'kepsilonphitf':
            self.foam_field_list = ['k',
                                'epsilon',
                                'T_t_ke',
                                'T_t_nut',
                                'T_k',
                                'U',
                                'gradp',
                                'gradk',
                                'gradv2',
                                'p',
                                'DUDt',
                                'wallDistance',
                                'nut',
                                'phit',
                                'f',
                                'S',
                                'Shat',
                                'R',
                                'Rhat',
                                'Av2',
                                'Ak',
                                'Av2hat',
                                'Akhat',
                                'gradU',
                                'skewness',
                                'C'
                                ]
            self.read_invariants_flag = True
            self.read_basis_tensors_flag = True
            self.read_q_flag = True
            self.read_lda_flag = True
            self.save_mesh_skewness = True

        if self.case_type == 'komegasst' or 'komega':
            self.foam_field_list = ['k',
                                'omega',
                                'epsilon',
                                'T_t_ke',
                                'T_t_nut',
                                'T_k',
                                'U',
                                'gradp',
                                'gradk',
                                'gradomega',
                                'p',
                                'DUDt',
                                'wallDistance',
                                'nut',
                                'S',
                                'Shat',
                                'R',
                                'Rhat',
                                'Ao',
                                'Ak',
                                'Aohat',
                                'Akhat',
                                'gradU',
                                'skewness',
                                'C'
                                ]
            self.read_invariants_flag = True
            self.read_basis_tensors_flag = True
            self.read_q_flag = True
            self.read_lda_flag = True
            self.save_mesh_skewness = True

        if self.case_type == 'LES':
            self.foam_field_list = [
                                'UMean',
                                'gradUMean',
                                'tauMean',
                                'SMean',
                                'RMean',
                                'kMean',
                                'kMean_tauMean',
                                'aMean',
                                'bMean',
                                'C'
                                ]
            self.read_invariants_flag = False
            self.read_basis_tensors_flag = False
            self.read_q_flag = False
            self.read_lda_flag = False
            self.save_mesh_skewness = False
            
        if self.case_type == 'DNS':
            self.foam_field_list = [
                                'U',
                                'gradU',
                                'TauDNS',
                                'S',
                                'R',
                                'k',
                                'a',
                                'b',
                                'C'
                                ]
            self.read_invariants_flag = False
            self.read_basis_tensors_flag = False
            self.read_q_flag = False
            self.read_lda_flag = False
            self.save_mesh_skewness = False

        self.assemble_full_foam_field_list()

    def assemble_full_foam_field_list(self):
        if self.read_invariants_flag:
            for i in range(47):
                self.foam_field_list.append(f'I1_{i+1}')
                self.foam_field_list.append(f'I2_{i+1}')
        if self.read_basis_tensors_flag:
            for i in range(10):
                self.foam_field_list.append(f'T{i+1}')
        if self.read_q_flag:
            for i in range(9):
                self.foam_field_list.append(f'q{i+1}')
        if self.read_lda_flag:
            for i in range(5):
                self.foam_field_list.append(f'lambda{i+1}')
        print('[dataFoam] Assembled full foam field list: ')
        print(self.foam_field_list)

    def writeFields(self):
        """Copies the foam case to the writeFields directory, changes startTime to latestTime, then calls the write_fields_application"""
        if not os.path.isdir(os.path.join(self.foam_parent_dir, 'writeFields')):
            os.mkdir(os.path.join(self.foam_parent_dir, 'writeFields'))

        self.writeFieldsDirectory = os.path.join(self.foam_parent_dir,'writeFields',self.case_name)
        if (self.write_fields_flag):
            print('[dataFoam] Writing new fields....')
            if (os.path.isdir(os.path.join(self.writeFieldsDirectory))):
                os.system(f'rm -rf {self.writeFieldsDirectory}')
            os.system('cp -rf '+self.directory+' '+self.writeFieldsDirectory)
            os.chdir(self.writeFieldsDirectory)
            changeFoamSystemDictEntry(os.path.join(self.writeFieldsDirectory),'controlDict','startFrom','latestTime')
            if self.save_mesh_skewness:
                print(f'[dataFoam] Running checkMesh....')
                os.system(f'checkMesh -writeFields skewness -time {self.endtime} > log.checkMesh')
            print('[dataFoam] Getting cell centres for the case....')
            self.C = get_cell_centres(self.writeFieldsDirectory)  
            print(f'[dataFoam] Running {self.writeFieldsApplication}....')
            os.system(f'{self.writeFieldsApplication} > log.writeFields')
        else:
            print('[dataFoam] Skipping writing fields....')
        return

    def saveDataset(self,dataset_prefix):
        """Read foam fields, save foam fields as numpy binaries"""
        self.foamdatatime = os.path.join(self.writeFieldsDirectory,str(self.endtime))
        self.dataset_prefix = dataset_prefix #+ '_'+self.case_name
        self.save_dir = os.path.join(self.data_save_path,dataset_prefix)

        print('[dataFoam] Reading fields.... ')
        for field in self.foam_field_list:
            foamfield = readFoamField(os.path.join(self.foamdatatime,field))
            if isinstance(foamfield,float):
                print(f'[dataFoam] Field {field} is all {foamfield}.')
                foamfield = np.ones(len(self.C))*foamfield
            print(f'[dataFoam] Saving {field}, with shape {foamfield.shape}')
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+field+'.npy'),foamfield)  
        return      