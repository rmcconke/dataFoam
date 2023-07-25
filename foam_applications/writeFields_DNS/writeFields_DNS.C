/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

This utility is a modified version of the simpleFoam solver, that does not solve any equations. 
It loads the OpenFOAM fields at the time given by startTime in the controlDict, and writes several fields as shown below.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    Info<< "Calculating extra DNS fields...." << nl << endl;

    dimensionedScalar k_small("k_small",dimensionSet(0,2,-2,0,0,0,0), 1E-10);
    k = 0.5*tr(tau); //TKE
    a = dev(tau); //Anisotropic part of Reynolds stress tensor
    b = a/(max(k_small,2.0*k)); //Non-dimensional anisotropy tensor

    // Velocity gradient and related tensors
    gradU = fvc::grad(U);
    gradU = gradU.T(); 	// Output the Jacobian, a more common form of the velocity gradient tensor
    S = symm(fvc::grad(U));
    R = -skew(fvc::grad(U)); // grad(UMean) produces the transpose of the Jacobian, RMean is defined based on the Jacobian, hence negative sign

    Info<< "Writing extra DNS fields...." << nl << endl;
    tau.write();
    a.write();
    b.write();
    k.write();
    gradU.write();
    S.write();
    R.write();
    Info<< "Finished writing extra DNS fields." << nl << endl;

    runTime.printExecutionTime(Info);
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
