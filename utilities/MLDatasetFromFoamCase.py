#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:22:26 2021

@author: ryley
"""

import os
import numpy as np
from scipy.interpolate import griddata
from utilities.foamIO.readFoam import readFoamTensor, readFoamSymmTensor, readFoamScalar, readFoamVector, get_cell_count, get_endtime, get_cell_centres, get_cell_volumes

class MLDatasetFromFoamCase: 
    # Class for a foam case that we want to store a numpy dataset for. 
    # The fields stored depend on the case_type.
    # RANS case types store the full set of invariants, while reference (LES/DNS)
    # only stores U/gradU/tau fields.
    
    def __init__(self,data_save_path,foam_parent_dir,case_name,case_type):
        # ML path - path to main ML folder
        # foam_parent_directory - e.g., the kepsilonphitf dataset directory
        # case_name - unique identifier
        # case_type: kepsilonphitf, komegasst, reference
        print('    Initializing MLDatasetFromFoamCase....')
        self.case_type = case_type
        self.foam_parent_dir = foam_parent_dir
        self.directory = os.path.join(self.foam_parent_dir,case_name)
        print('        Case directory: '+self.directory)

        self.case_name = case_name
        self.data_save_path = data_save_path
        self.cells = get_cell_count(self.directory)
        self.endtime = get_endtime(self.directory)
        print('    Getting cell centres for the case....')
        self.Cx, self.Cy, self.Cz = get_cell_centres(self.directory)        
        
        if self.case_type == 'kepsilonphitf':
            self.foam_scalars_list = ['k',
                                'epsilon',
                                'T_t',
                                'T_k',
                                'Cx',
                                'Cy',
                                'Cz',
                                'Ux',
                                'Uy',
                                'Uz',
                                'gradpx',
                                'gradpy',
                                'gradpz',
                                'gradkx',
                                'gradky',
                                'gradkz',
                                'gradv2x',
                                'gradv2y',
                                'gradv2z',
                                'p',
                                'DUDtx',
                                'DUDty',
                                'DUDtz',
                                'wallDistance',
                                'phit',
                                'f',
                                ]
                
            self.foam_symmtensors_list = ['S',
                                     'Shat',
                                     ]
                
            self.foam_tensors_list = ['R',
                                'Rhat',
                                'Av2',
                                #'Ap',
                                'Ak',
                                'Av2hat',
                                #'Aphat',
                                'Akhat',
                                'gradU',
                                ]
            self.read_invariants_flag = True

        if self.case_type == 'komegasst':
            self.foam_scalars_list = ['k',
                                'epsilon',
                                'T_t',
                                'T_k',
                                'Cx',
                                'Cy',
                                'Cz',
                                'Ux',
                                'Uy',
                                'Uz',
                                'gradpx',
                                'gradpy',
                                'gradpz',
                                'gradkx',
                                'gradky',
                                'gradkz',
                                'omega',
                                'gradomegax',
                                'gradomegay',
                                'gradomegaz',
                                'p',
                                'DUDtx',
                                'DUDty',
                                'DUDtz',
                                'wallDistance',
                                ]
                
            self.foam_symmtensors_list = ['S',
                                     'Shat',
                                     ]
                
            self.foam_tensors_list = ['R',
                                'Rhat',
                                'Ap',
                                'Ak',
                                'Aphat',
                                'Akhat',
                                'gradU',
                                ]   
            self.read_invariants_flag = True

            
        if self.case_type == 'reference':
            self.foam_scalars_list = [
                                'Cx',
                                'Cy',
                                'Cz',
                                'Ux',
                                'Uy',
                                'Uz',
                                ]
                
            self.foam_symmtensors_list = ['tau']
            self.foam_tensors_list = ['gradU'] 
            self.read_invariants_flag = False
            
    def writeFields(self):
        self.writeFieldsDirectory = os.path.join(self.foam_parent_dir,'02-writeFields',self.case_name)
        if not os.path.isdir(self.writeFieldsDirectory):
            os.system('cp -rf '+self.directory+' '+self.writeFieldsDirectory)
        
        if not os.path.isfile(os.path.join(self.writeFieldsDirectory,'Write.complete')):
            os.system('cp -r '+os.path.join(self.foam_parent_dir,'02-writeFields','Runwrite')+' '+self.writeFieldsDirectory)
            os.chdir(self.writeFieldsDirectory)
            os.system('./Runwrite')
            os.system('touch Write.complete')
        else:
            print('    Looks like fields have already been written, moving on....')
            
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
            
    def saveDataset(self,dataset_prefix):
        self.foamdatatime = os.path.join(self.writeFieldsDirectory,str(self.endtime))
        self.foam_scalars_dict = dict.fromkeys(self.foam_scalars_list)
        self.foam_tensors_dict = dict.fromkeys(self.foam_tensors_list)
        self.foam_symmtensors_dict = dict.fromkeys(self.foam_symmtensors_list)
        
        self.dataset_prefix = dataset_prefix #+ '_'+self.case_name
        self.save_dir = os.path.join(self.data_save_path,dataset_prefix)


        print('Reading scalars.... ')
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
    
    def get_cell_volumes(self):
        self.CV = get_cell_volumes(self.directory)
 
