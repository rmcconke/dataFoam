#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:25:35 2021

@author: ryley
"""
import numpy as np

def writeFoam_U_DUCT(filename, U):
    # Writes the aperp file
    cells = len(U)
    leftbracket = np.repeat('(',cells)
    rightbracket = np.repeat(')',cells)
    array_write = np.column_stack((leftbracket,U,rightbracket))
    print('[dataFoam] Writing U to file '+filename)
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

internalField nonuniform List<vector>
""")
    file.write(str(cells) + """ 
(
""")
    np.savetxt(file, array_write,fmt='%s')
    file.write(""");
boundaryField
{
    inlet
    {
        type            cyclic;
    }

    outlet
    {
        type            cyclic;
    }

    walls
    {
        type            noSlip;
    }
}
// ************************************************************************* //""")
    file.close()

def writeFoam_TauDNS_DUCT(filename, tau):
    # Writes the tau file, assumes tau is column stacked
    cells = len(tau)
    leftbracket = np.repeat('(',cells)
    rightbracket = np.repeat(')',cells)
    array_write = np.column_stack((leftbracket,tau,rightbracket))
    print('[dataFoam] Writing TauDNS to file '+filename)
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
    class       volSymmTensorField;
    object      TauDNS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField nonuniform List<symmTensor>
""")
    file.write(str(cells) + """
(
""")
    np.savetxt(file, array_write,fmt='%s')
    file.write(""");
boundaryField
{
    inlet
    {
        type            cyclic;
    }

    outlet
    {
        type            cyclic;
    }

    ".*"
    {
        type            fixedValue;
        value           uniform (0 0 0 0 0 0);
    }
}
    // ************************************************************************* //""")
    file.close()
