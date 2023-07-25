#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 17:25:35 2021

@author: ryley
"""
import numpy as np
from numpy import inf

import os

def foam_tensor(field):
    foam_shaped = np.column_stack((field[:,0,0],field[:,0,1],field[:,0,2],field[:,1,0],field[:,1,1],field[:,1,2],field[:,2,0],field[:,2,1],field[:,2,2]))
    return foam_shaped

def writeFoam_anyfield_CUBE(field,filename):
    # Writes a generic scalar field with zero dimensions, for visualization
    field_name=os.path.basename(filename)
    #field[field == -inf] = -1E6
    #field[field == inf] = 1E6
    cells = len(field)
    print('[dataFoam] Writing {} to file {}'.format(field_name,filename))
    nan_count = np.count_nonzero(np.isnan(field))
    if nan_count > 0:
        print('\n[dataFoam] WARNING! Found '+str(nan_count)+' NaN values in the field. Replacing these with 1E10 when writing to foam, but will still be present in saved array.')
        field=np.nan_to_num(field, copy=True, nan=1E10)

    fieldclass = 'volScalarField'
    listtype='scalar'
    zero_value = '0'
    if field.ndim > 1:
        if field.ndim == 2: 
            fieldclass = 'volVectorField'
            listtype='vector'
            zero_value = '(0 0 0)'

        elif field.ndim == 3:
            field = foam_tensor(field)
            fieldclass = 'volTensorField'
            listtype='tensor'
            zero_value = '(0 0 0 0 0 0 0 0 0)'

        leftbracket = np.repeat('(',len(field))
        rightbracket = np.repeat(')',len(field))
        field = np.column_stack((leftbracket,field,rightbracket))

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
    class       """+fieldclass+""";
    object      """+field_name+""";
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField nonuniform List<"""+listtype+""">
""")
    file.write(str(cells) + """ (
        """)
    np.savetxt(file, field,fmt='%s')
    file.write(""");
boundaryField
{
    INLET
    {
        type            calculated;
        value           uniform """+zero_value+""";
    }
    CUBE1
    {
        type            fixedValue;
        value           uniform """+zero_value+""";
    }
    CUBE2
    {
        type            fixedValue;
        value           uniform """+zero_value+""";
    }
    WALL
    {
        type            fixedValue;
        value           uniform """+zero_value+""";
    }           
    RIGHT
    {
        type            cyclic;
    }
    LEFT
    {
        type            cyclic;
    }
    TOP
    {
        type            zeroGradient;
    }
    OUTLET
    {
        type            zeroGradient;
    }


}
// ************************************************************************* //""")
    file.close()

def writeFoam_ap_CUBE(filename, aperp):
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
        INLET
        {
            type            calculated;
            value           uniform (0 0 0 0 0 0);
        }
        CUBE1
        {
            type            fixedValue;
            value           uniform (0 0 0 0 0 0);
        }
        CUBE2
        {
            type            fixedValue;
            value           uniform (0 0 0 0 0 0);
        }
        WALL
        {
            type            fixedValue;
            value           uniform (0 0 0 0 0 0);
        }           
        RIGHT
        {
            type            cyclic;
        }
        LEFT
        {
            type            cyclic;
        }
        TOP
        {
            type            zeroGradient;
        }
        OUTLET
        {
            type            zeroGradient;
        }

}
// ************************************************************************* //""")
    file.close()

def writeFoam_nut_L_CUBE(filename,nut_L):
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
    INLET
    {
        type            calculated;
        value           uniform 0;
    }
    CUBE1
    {
        type            fixedValue;
        value           uniform 0;
    }
    CUBE2
    {
        type            fixedValue;
        value           uniform 0;
    }
    WALL
    {
        type            fixedValue;
        value           uniform 0;
    }           
    RIGHT
    {
        type            cyclic;
    }
    LEFT
    {
        type            cyclic;
    }
    TOP
    {
        type            zeroGradient;
    }
    OUTLET
    {
        type            zeroGradient;
    }


}
// ************************************************************************* //""")
    file.close()
    
def writeFoam_pfv_CUBE(filename, pfv):
    pfv[pfv == -inf] = -1E6
    pfv[pfv == inf] = 1E6

    # Writes the pfv file
    cells = len(pfv)
    #aperp = np.column_stack((aperp[:,0,0],aperp[:,0,1],aperp[:,0,2],
    #                                   aperp[:,1,1],aperp[:,1,2],
    #                                                -aperp[:,0,0]-aperp[:,1,1]))
    leftbracket = np.repeat('(',len(pfv))
    rightbracket = np.repeat(')',len(pfv))
    array_write = np.column_stack((leftbracket,pfv,rightbracket))
    print('[dataFoam] Writing pfv to file '+filename)
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
        object      pfv;
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    dimensions      [0 1 -2 0 0 0 0];
    
    internalField nonuniform List<vector>
    """)
    file.write(str(cells) + """ (
        """)
    np.savetxt(file, array_write,fmt='%s')
    file.write(""");
    boundaryField
    {
            INLET
            {
                type            calculated;
                value           uniform (0 0 0);
            }
            CUBE1
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            CUBE2
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            WALL
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }           
            RIGHT
            {
                type            cyclic;
            }
            LEFT
            {
                type            cyclic;
            }
            TOP
            {
                type            zeroGradient;
            }
            OUTLET
            {
                type            zeroGradient;
            }
    
    }
    // ************************************************************************* //""")
    file.close()