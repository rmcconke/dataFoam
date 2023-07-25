#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:25:35 2021

@author: ryley
"""
import numpy as np
import os
import re
import Ofpp

def readFoamField(file):
    field = Ofpp.parse_internal_field(file)
    if isinstance(field, float):
        return field
    elif field.ndim > 1: 
        if field.shape[1] == 6:
            field = reshape_symmTensor(field)
        elif field.shape[1] == 9:
            field = reshape_tensor(field)
    return field

def reshape_symmTensor(Tensor):
    Tensor = np.stack((Tensor[:,0], Tensor[:,1], Tensor[:,2], Tensor[:,1], Tensor[:,3], Tensor[:,4], Tensor[:,2], Tensor[:,4],Tensor[:,5]),axis=1)
    Tensor = reshape_tensor(Tensor)
    return Tensor

def reshape_tensor(Tensor):
    Tensor = Tensor.reshape((len(Tensor),3,3))
    return Tensor

def get_cell_count(foam_directory):
    # Returns scalar number of cells for the case at foam_directory.
    # Uses a trick to read a specific part of the polyMesh/neighbour file.
    try: 
        with open(os.path.join(foam_directory,'constant','polyMesh','neighbour')) as file:
            lines = file.readlines()
            for line in lines:
                if re.search(r'\s* note',line):
                    cells = int(re.search(r'nCells:(.+?) ',line).group(1))
                    break
        print('[dataFoam] Found a mesh with number of cells: '+str(cells))
    except:
        raise LookupError('[dataFoam] Could not find a mesh for the foam case '+foam_directory)
    return cells

def get_endtime(foam_directory):
    # Returns the largest reconstructed time step as an int for a foam_directory.
    try: 
        folders_list = next(os.walk(foam_directory))[1]
        endtime=0
        folders_list = [ x for x in folders_list if "processor" not in x ]
        for folder in folders_list:
            x = re.findall("\\d", folder)
            result = "".join(x)
            if result != '':
                if int(result) > endtime:
                    endtime=int(result)
                
        print('[dataFoam] Found endtime: '+ str(endtime))
    except: 
        raise LookupError('[dataFoam] Could not get endtime for the foam case '+foam_directory)
    return endtime

def get_cell_centres(foam_directory):
    # Returns three vectors, Cx, Cy, Cz. Runs writeCellCentres then reads the files.
    logfile = os.path.join(foam_directory,'log.writeCellCentres')
    print(f'[dataFoam] Running writeCellCentres....')
    os.system(f'postProcess -case {foam_directory} -func writeCellCentres > {logfile}')
    cell_centres = readFoamField(os.path.join(foam_directory,'0/C'))
    #Cx, Cy, Cz = cell_centres[:,0], cell_centres[:,1], cell_centres[:,2]  
    return cell_centres
    
def readLineofFile(file,line):
    # Returns a line of a file.
    a_file = open(file)
    lines_to_read = [line]    
    for position, line in enumerate(a_file):        
        if position in lines_to_read:
            #print('Reading line: ')
            #print(line)
            read = line
    a_file.close()
    print('[dataFoam] Read line '+str(line) + ' of file '+ file + ' : '+ str(read))
    return read
