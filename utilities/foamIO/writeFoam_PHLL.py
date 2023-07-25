#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:25:35 2021

@author: ryley
"""
import numpy as np

def writeFoam_genericscalar_PHLL(filename,field_name,field):
    # Writes a generic scalar field with zero dimensions, for visualization
    cells = len(field)
    print('[dataFoam] Writing {} to file {}'.format(field_name,filename))
    nan_count = np.count_nonzero(np.isnan(field))
    if nan_count > 0:
        print('\n[dataFoam] WARNING! Found '+str(nan_count)+' NaN values in the field. Replacing these with 1E10 when writing to foam, but will still be present in saved array.')
        field=np.nan_to_num(field, copy=True, nan=1E10)
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
    object      """+field_name+""";
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField nonuniform List<scalar>
""")
    file.write(str(cells) + """ (
""")
    np.savetxt(file, field)
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
        value           uniform 0;
    }    
// ************************************************************************* //""")
    file.close()

def writeFoam_genericsymmTensor_PHLL(filename, field_name, tensor):
    # Writes the aperp file
    cells = len(tensor)
    tensor = np.column_stack((tensor[:,0,0],tensor[:,0,1],tensor[:,0,2],
                                       tensor[:,1,1],tensor[:,1,2],
                                                    tensor[:,2,2]))
    leftbracket = np.repeat('(',len(tensor))
    rightbracket = np.repeat(')',len(tensor))
    array_write = np.column_stack((leftbracket,tensor,rightbracket))
    print(f'[dataFoam] Writing {field_name} to file {filename}')
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
    object      """+field_name+""";
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField nonuniform List<symmTensor>
""")
    file.write(str(cells) + """ (
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


def writeFoam_ap_PHLL(filename, aperp):
    # Writes the aperp file
    cells = len(aperp)
    aperp = np.column_stack((aperp[:,0,0],aperp[:,0,1],aperp[:,0,2],
                                       aperp[:,1,1],aperp[:,1,2],
                                                    -aperp[:,0,0]-aperp[:,1,1]))
    leftbracket = np.repeat('(',len(aperp))
    rightbracket = np.repeat(')',len(aperp))
    array_write = np.column_stack((leftbracket,aperp,rightbracket))
    print('[dataFoam] Writing aperp to file '+filename)
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
    object      aperp;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField nonuniform List<symmTensor>
""")
    file.write(str(cells) + """ (
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

def writeFoam_nut_L_PHLL(filename,nut_L):
    # Writes the nut_L file
    cells = len(nut_L)
    print('[dataFoam] Writing nut_L to file {}'.format(filename))
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
    object      nut_L;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField nonuniform List<scalar>
""")
    file.write(str(cells) + """ (
""")
    np.savetxt(file, nut_L, fmt='%s')
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
        value           uniform 0;
    }
}
// ************************************************************************* //""")
    file.close()
