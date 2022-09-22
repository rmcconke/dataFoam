#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 14:28:11 2021

@author: ryley
"""
import fileinput
import sys
import os

def changeFoamSystemDictEntry(foam_dir,dict_name,entry,value):
    # Changes an entry in the system/dict_name to be value
    print('Setting new foam '+entry+ ' :' + str(value)+ ' in directory '+foam_dir +' '+dict_name)
    for line in fileinput.input(os.path.join(foam_dir,'system',dict_name), inplace=True):
        if line.strip().startswith(entry):
            line = entry+'         '+ str(value) + ';\n'
        sys.stdout.write(line)
    return
