#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:25:35 2021

@author: ryley
"""
import numpy as np
import os
import re

def readFoamVector(file, cells=None):
    # Returns a numpy vector of shape (cells,3) for the internal vector field
    # at the input file. If something goes wrong, returns a uniform zero vector.
    print('    Reading vector: '+file)
    cleanFoamFile(file)
    if cells == None:
        cells= get_cell_count(os.path.dirname(os.path.dirname(file)))
    try:
        Vector = np.genfromtxt(file+'_cleaned',skip_header=23,max_rows=cells)
        removeCleanedFoamFile(file)
    except: 
        print('--->ERROR! Could not read vector: '+file+', setting to uniform zero')
        Vector = np.zeros(cells*3)
    Vector = Vector.reshape((cells,3))
    return Vector

def readFoamTensor(file,cells=None):
    # Returns a numpy tensor of shape (cells,3,3) for the internal tensor field
    # at the input file. If something goes wrong, returns a uniform zero tensor.
    print('    Reading tensor: '+file)
    cleanFoamFile(file)
    if cells == None:
        cells= get_cell_count(os.path.dirname(os.path.dirname(file)))
    try:
        Tensor = np.genfromtxt(file+'_cleaned',skip_header=23,max_rows=cells)
        removeCleanedFoamFile(file)
    except: 
        print('--->ERROR! Could not read tensor: '+file+', setting to uniform zero')
        Tensor = np.zeros(cells*9)
    Tensor = Tensor.reshape((cells,3,3))
    return Tensor

def readFoamSymmTensor(file,cells=None):
    # Returns a numpy tensor of shape (cells,3,3) for the internal symmTensor field
    # at the input file. If something goes wrong, returns a uniform zero symmTensor.
    print('    Reading symmTensor: '+file)
    cleanFoamFile(file)
    if cells == None:
        cells= get_cell_count(os.path.dirname(os.path.dirname(file)))
    try:
        Tensor = np.genfromtxt(file+'_cleaned',skip_header=23,max_rows=cells)
        removeCleanedFoamFile(file)
    except: 
        print('--->ERROR! Could not read symmTensor: '+file+', setting to uniform zero')
        Tensor = np.zeros(cells*6)
    Tensor = np.stack((Tensor[:,0], Tensor[:,1], Tensor[:,2], Tensor[:,1], Tensor[:,3], Tensor[:,4], Tensor[:,2], Tensor[:,4],Tensor[:,5]),axis=1)
    Tensor = Tensor.reshape((cells,3,3))
    return Tensor

def readFoamScalar(file,cells=None):
    # Returns a numpy scalar of shape (cells) for the internal scalar field
    # at the input file. If something goes wrong, returns a uniform zero scalar.
    print('    Reading scalar: '+file)
    if cells == None:
        cells= get_cell_count(os.path.dirname(os.path.dirname(file)))
    try:
        scalar=np.genfromtxt(file,skip_header=23,max_rows=cells)
    except: 
        print('--->ERROR! Could not read scalar: '+file+', setting to uniform zero')
        scalar = np.zeros(cells)
    return scalar

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
        print('    Found a mesh with number of cells: '+str(cells))
    except:
        raise LookupError('Could not find a mesh for the foam case '+foam_directory)
    return cells

def get_endtime(foam_directory):
    # Returns the largest reconstructed time step as an int for a foam_directory.
    try: 
        folders_list = next(os.walk(foam_directory))[1]
        endtime=0
        for folder in folders_list:
            x = re.findall("\\d", folder)
            result = "".join(x)
            if result != '':
                if int(result) > endtime:
                    endtime=int(result)
                
        print('    Found endtime: '+ str(endtime))
    except: 
        raise LookupError('Could not get endtime for the foam case '+foam_directory)
    return endtime

def get_cell_centres(foam_directory):
    # Returns three vectors, Cx, Cy, Cz. Runs writeCellCentres then deletes
    # the Cx, Cy, Cz files.
    try:
        os.system('postProcess -case ' + foam_directory 
                  + ' -func writeCellCentres >' + foam_directory+'/log.writeCellCentres')
        Cx = readFoamScalar(os.path.join(foam_directory,str(0),'Cx'))
        Cy = readFoamScalar(os.path.join(foam_directory,str(0),'Cy'))
        Cz = readFoamScalar(os.path.join(foam_directory,str(0),'Cz'))
        os.system('rm -r '+foam_directory+'/0/Cx')
        os.system('rm -r '+foam_directory+'/0/Cy')
        os.system('rm -r '+foam_directory+'/0/Cz')
        os.system('rm -r '+foam_directory+'/0/C')
    except:
        raise LookupError('Could not get cell centres for the foam case '+foam_directory)
    return Cx, Cy, Cz

def get_cell_volumes(foam_directory):
    # Returns the cell volumes for foam_directory as a vector. 
    # Runs writeCellCentres then deletes the V file.
    try: 
        os.system('postProcess -func writeCellVolumes > log.writeCellVolumes')
        V = readFoamScalar(os.path.join(foam_directory,str(0),'V'))
        os.system('rm -r '+foam_directory+'/0/V')
    except:
        raise LookupError('Could not get cell volumes for the foam case '+foam_directory)
    return V
    
def cleanFoamFile(file):
    # Uses sed to strip the leading/trailing brackets from a file, creating
    # a new file with _cleaned appended. Used for reading vectors/tensors
    os.system("sed 's/(//g' " + file +' > ' + file+'_cleaned')
    os.system("sed -i 's/)//g' " + file+'_cleaned')
    
def removeCleanedFoamFile(file):
    os.system('rm -r '+file+'_cleaned')

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
    print('    Read line '+str(line) + ' of file '+ file + ' : '+ str(read))
    return read
