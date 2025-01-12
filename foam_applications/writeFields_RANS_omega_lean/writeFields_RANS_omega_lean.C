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

    Info<< "Calculating extra RANS fields...." << nl << endl;

	epsilon = 0.09*k*omega;
	T_t_ke = k/epsilon; //Time scale based on k and epsilon

	gradU = fvc::grad(U);
	gradU = gradU.T(); 	// Output the Jacobian, a more common form of the velocity gradient tensor

	turbR = (2.0/3.0)*I*k - (nut)*dev(twoSymm(fvc::grad(U)));
	divturbR=fvc::div(turbR);

	DUDt = U & gradU.T(); 	// Checked with component wise expression for the convective derivative

	S = symm(fvc::grad(U));
	R = -skew(fvc::grad(U)); // grad(U) produces the transpose of the Jacobian, R is defined based on the Jacobian, hence negative sign

	// Pressure gradient
	gradp = fvc::grad(p);
	gradomega = fvc::grad(omega);
	gradepsilon = fvc::grad(epsilon);
	gradnut = fvc::grad(nut);

	// TKE gradient
	gradk = fvc::grad(k);

	// Other turbulence fields
	wallDistance = wallDist(mesh).y();

	// Some q's from ref Kaandorp 2020, 10.1016/j.compfluid.2020.104497
	q1 = 0.5 * (mag(R)*mag(R) - mag(S)*mag(S)) /(max(mag(S)*mag(S),SMALL_S2)); // Ratio of excess rotation rate to strain rate
	q2 = min((sqrt(k)*wallDistance/(50.0*turbulence->nu())),2.0);  // Wall distance based Reynolds Number
	q3 = T_t_ke * mag(S); // Ratio of turbulent time scale to mean strain time scale
	q4 = mag(turbR)/(k); //Ratio of total to 1/2 * normal Reynolds stresses (TKE)

	// Other q's (i.e., heuristic input features) can also be added to this code, by adding a new field in createFields.H, and assigning values here

    Info<< "Writing extra RANS fields...." << nl << endl;

	// Write fields
	turbR.write();
	divturbR.write();
	q1.write();
	q2.write();
	q3.write();
	q4.write();
	S.write();
	R.write();
	gradp.write();
	gradk.write();
	gradU.write();
	gradepsilon.write();
	gradnut.write();
	gradomega.write();
	p.write();
	epsilon.write();
	runTime.write();
	DUDt.write();
	wallDistance.write();
	runTime.printExecutionTime(Info);
    
	Info<< "Finished writing extra RANS fields." << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
