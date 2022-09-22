#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:55:19 2021

@author: ryley
"""

def writeFoamktheta_CBFS(filename,cells,k):
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
    class       volScalarField;
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
dimensions      [0 2 -2 0 0 0 0];
    
internalField   nonuniform List<scalar>
     """)
    file.write(str(cells) + """ (
        
        """)
    for i in range(cells):
        file.write(str(k[i]) +'\n')   
    file.write(""");
    boundaryField
    {
     "(top|bottom)"
    {
        type            uniform 0;
    }
    "(inlet|front|back|outlet)"
    {
        type            zeroGradient;
    }

    }          

     
    // ************************************************************************* //""")
    file.close()