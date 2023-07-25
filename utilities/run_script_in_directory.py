#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

def run_script_in_directory(script,directory):
    subprocess.check_call(['cp', '-r', script, directory])
    subprocess.check_call('sh '+directory+'/'+script, shell=True)
    return