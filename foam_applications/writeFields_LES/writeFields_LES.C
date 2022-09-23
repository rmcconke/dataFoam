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
    /*argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );
	*/
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    //#include "initContinuityErrs.H"

    //turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

		Info<< "Calculating fields for LES" << nl << endl;
        // Reynolds stress tensor and related tensors
		//tauMean = UPrime2Mean+turbulence->R();
        // Note: currently assumes eddy viscosity SGS
        subgrid_tauMean = ((2.0/3.0)*I)*kMean_model - (nutMean)*dev(twoSymm(fvc::grad(UMean))),
		tauMean = UPrime2Mean+subgrid_tauMean;
        aMean = dev(tauMean);
        bMean = aMean/(tr(tauMean));
        kMean = 0.5*tr(tauMean);

        // Velocity gradient and related tensors
		gradUMean = fvc::grad(UMean);
		gradUMean = gradUMean.T(); 	// Output the Jacobian, a more common form of the velocity gradient tensor
       	SMean = symm(fvc::grad(UMean));
		RMean = -skew(fvc::grad(UMean)); // grad(UMean) produces the transpose of the Jacobian, RMean is defined based on the Jacobian, hence negative sign

		Info<< "Writing fields for LES...." << nl << endl;
        subgrid_tauMean.write();
		tauMean.write();
        aMean.write();
        bMean.write();
        kMean.write();
        gradUMean.write();
        SMean.write();
        RMean.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
