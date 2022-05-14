#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:55:19 2021

@author: ryley
"""

def writeFoamUtheta_CNDV(filename,cells,Ux,Uy,Uz):
    #filename = directory+'/0/Utheta'
    print('Writing interpolated Utheta' + ' to file '+filename)
    file = open(filename,'w')
    file.write("""/*--------------------------------*- C++ -*----------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  v2006                                 |
    |   \\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       volVectorField;
        object      U;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    dimensions      [0 1 -1 0 0 0 0];
    
    internalField   nonuniform List<vector>
     """)
    file.write(str(cells) + """ (
        """)
    for i in range(cells):
        file.write('('+str(Ux[i])+' '+str(Uy[i])+ ' '+str(Uz[i]) +')\n')   
    file.write(""");
    boundaryField
    {
    "(bottom|top)"
    {
        type            noSlip;
    }
    outlet
    {
        type            zeroGradient;
    }
    inlet
    {
        type            zeroGradient;
    }

    }          

     
    // ************************************************************************* //""")
    file.close()