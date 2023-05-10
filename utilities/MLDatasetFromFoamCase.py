#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

import os
import numpy as np
from dataFoam.utilities.foamIO.readFoam import readFoamField, get_endtime, get_cell_centres

class MLDatasetFromFoamCase: 
    # Class for a foam case that we want to store a numpy dataset for. 
    # The fields stored depend on the case_type.
    # RANS case types store the full set of invariants, while reference (LES/DNS)
    # only stores U/gradU/tau fields.
    
    def __init__(self,data_save_path,foam_parent_dir,case_name,case_type,overwrite_flag=True):
        # ML path - path to main ML folder
        # foam_parent_directory - e.g., the kepsilonphitf dataset directory
        # case_name - unique identifier
        # case_type: kepsilonphitf, komegasst, reference
        print('[dataFoam] Initializing MLDatasetFromFoamCase....')
        self.case_type = case_type
        self.foam_parent_dir = foam_parent_dir
        self.directory = os.path.join(self.foam_parent_dir,case_name)
        print('[dataFoam] Case directory: '+self.directory)

        self.case_name = case_name
        self.data_save_path = data_save_path
        self.overwrite_flag=overwrite_flag

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

        if self.case_type == 'komegasst':
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
                                'Ap',
                                'Ak',
                                'Aphat',
                                'Akhat',
                                'gradU',
                                'skewness',
                                'C'
                                ]
            self.writeFieldsApplication='writeFields_komegasstMean'
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
            self.writeFieldsApplication='writeFields_LES'
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
            for i in range(4):
                self.foam_field_list.append(f'q{i+1}')
        if self.read_lda_flag:
            for i in range(5):
                self.foam_field_list.append(f'lambda{i+1}')
        print('[dataFoam] Assembled full foam field list: ')
        print(self.foam_field_list)

    def writeFields(self):
        if not os.path.isdir(os.path.join(self.foam_parent_dir, '02-writeFields')):
            os.mkdir(os.path.join(self.foam_parent_dir, '02-writeFields'))

        self.writeFieldsDirectory = os.path.join(self.foam_parent_dir,'02-writeFields',self.case_name)
        if (self.overwrite_flag):
            print('[dataFoam] Writing new fields....')
            if (os.path.isdir(os.path.join(self.writeFieldsDirectory))):
                os.system(f'rm -rf {self.writeFieldsDirectory}')
            os.system('cp -rf '+self.directory+' '+self.writeFieldsDirectory)
            #os.system('cp -r '+os.path.join(self.foam_parent_dir,'02-writeFields','Runwrite')+' '+self.writeFieldsDirectory)
            os.chdir(self.writeFieldsDirectory)
            if self.save_mesh_skewness:
                os.system(f'checkMesh -writeFields skewness -time {self.endtime} | tee log.checkMesh')
            print('[dataFoam] Getting cell centres for the case....')
            self.C = get_cell_centres(self.writeFieldsDirectory)     
            os.system(f'{self.writeFieldsApplication} | tee log.writeFields')
        else:
            print('[dataFoam] Skipping writing fields....')

    def saveDataset(self,dataset_prefix):
        self.foamdatatime = os.path.join(self.writeFieldsDirectory,str(self.endtime))
        self.dataset_prefix = dataset_prefix #+ '_'+self.case_name
        self.save_dir = os.path.join(self.data_save_path,dataset_prefix)

        print('[dataFoam] Reading fields.... ')
        for field in self.foam_field_list:
            foamfield = readFoamField(os.path.join(self.foamdatatime,field))
            print(f'[dataFoam] Saving {field}, with shape {foamfield.shape}')
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+field+'.npy'),foamfield)        

"""
        for foam_scalar in self.foam_scalars_dict:
            self.foam_scalars_dict[foam_scalar] = readFoamScalar(os.path.join(self.foamdatatime,foam_scalar),cells=self.cells)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+foam_scalar+'.npy'),self.foam_scalars_dict[foam_scalar])
            
        print('Reading tensors....\n ')
        for foam_tensor in self.foam_tensors_dict:
            self.foam_tensors_dict[foam_tensor] = readFoamTensor(os.path.join(self.foamdatatime,foam_tensor),cells=self.cells)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+foam_tensor+'.npy'),self.foam_tensors_dict[foam_tensor]) 

        print('Reading symmtensors....\n ')

        for foam_symmtensor in self.foam_symmtensors_dict:
            self.foam_symmtensors_dict[foam_symmtensor] = readFoamSymmTensor(os.path.join(self.foamdatatime,foam_symmtensor),cells=self.cells)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+foam_symmtensor+'.npy'),self.foam_symmtensors_dict[foam_symmtensor]) 

        if self.read_invariants_flag:
            print('Reading invariants....\n ')
            self.read_invariants()
            print('Reading basis tensors....\n ')
            self.read_basis_tensors()
            print('Reading q....\n ')
            self.read_q()
            print('Reading lambdas....\n ')
            self.read_lda()
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_q.npy'), self.q)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_lambda.npy'), self.lda)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_Tensors.npy'), self.basis_tensors)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_I1.npy'), self.I1)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_I2.npy'), self.I2)
            
        if self.case_type == 'reference_mapped':
            pfv_x = self.foam_scalars_dict['gradaperpxx_x'] + self.foam_scalars_dict['gradaperpxy_x'] + self.foam_scalars_dict['gradaperpxz_x']
            pfv_y = self.foam_scalars_dict['gradaperpxy_y'] + self.foam_scalars_dict['gradaperpyy_y'] + self.foam_scalars_dict['gradaperpyz_y']
            pfv_z = self.foam_scalars_dict['gradaperpxz_z'] + self.foam_scalars_dict['gradaperpyz_z'] + self.foam_scalars_dict['gradaperpzz_z']
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+'pfv_x'+'.npy'),pfv_x)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+'pfv_y'+'.npy'),pfv_y)
            np.save(os.path.join(self.data_save_path,self.dataset_prefix+'_'+'pfv_z'+'.npy'),pfv_z)

    
    def get_cell_volumes(self):
        self.CV = get_cell_volumes(self.directory)
 
    def read_basis_tensors(self):
        bigT = np.empty((self.cells,10,3,3))
        for TN in range(0,10):
            tensorname = 'T' +str(TN+1)
            if TN == 0:
                bigT[:,TN,:,:] = readFoamSymmTensor(os.path.join(self.foamdatatime,tensorname),cells=self.cells)[:,:,:]
            else:
                bigT[:,TN,:,:] = readFoamTensor(os.path.join(self.foamdatatime,tensorname),cells=self.cells)[:,:,:]
        self.basis_tensors = bigT
    
    def read_invariants(self):
        self.I1 = np.empty((self.cells,47))
        for i in range(47):
            self.I1[:,i] = readFoamScalar(os.path.join(self.foamdatatime,'I1_'+str(i+1)),cells=self.cells)
        self.I2 = np.empty((self.cells,47))
        for i in range(47):
            self.I2[:,i] = readFoamScalar(os.path.join(self.foamdatatime,'I2_'+str(i+1)),cells=self.cells)
    
    def read_q(self):
        self.q = np.empty((self.cells,4))
        for i in range(4):
            self.q[:,i] = readFoamScalar(os.path.join(self.foamdatatime,'q'+str(i+1)),cells=self.cells)
            
    def read_lda(self):
        self.lda = np.empty((self.cells,5))
        for i in range(5):
            self.lda[:,i] = readFoamScalar(os.path.join(self.foamdatatime,'lambda'+str(i+1)),cells=self.cells)    
""" 