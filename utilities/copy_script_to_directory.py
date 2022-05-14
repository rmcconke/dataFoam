#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 16:30:50 2021

@author: ryley
"""
import subprocess

def copy_script_to_directory(script,directory):
    subprocess.check_call(['cp', '-r', script, directory])
    return
