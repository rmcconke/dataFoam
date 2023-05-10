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

		Info<< "Calculating fields for ML" << nl << endl;

		// Time scales
		epsilon = 0.09*k*omega;
		T_t_ke = k/epsilon;
		T_t_nut = nut/k;
		T_k = sqrt(turbulence->nu()/epsilon);
		C_Ak = nut/(sqrt(k)*k);
		
		// Velocity derivatives
		gradU = fvc::grad(U);
		gradU = gradU.T(); 	// Output the Jacobian, a more common form of the velocity gradient tensor
		DUDt = U & gradU.T(); 	// Checked with component wise expression for the convective derivative

       	S = symm(fvc::grad(U));
		R = -skew(fvc::grad(U)); // grad(U) produces the transpose of the Jacobian, R is defined based on the Jacobian, hence negative sign

		Shat = T_t_nut * S;
		Rhat = T_t_nut * R;

		// Pressure gradient
		gradp = fvc::grad(p);
		Ap.replace(1, -gradp.component(2));
		Ap.replace(2,  gradp.component(1));
		Ap.replace(3,  gradp.component(2));
		Ap.replace(5, -gradp.component(0));
		Ap.replace(6, -gradp.component(1));
		Ap.replace(7,  gradp.component(0));
		//Old normalization: Aphat = Ap/(max(mag(DUDt),SMALL_CONVDER));
		Aphat = C_Ak*Ap;

		// TKE gradient
		gradk = fvc::grad(k);
		Ak.replace(1, -gradk.component(2));
		Ak.replace(2,  gradk.component(1));
		Ak.replace(3,  gradk.component(2));
		Ak.replace(5, -gradk.component(0));
		Ak.replace(6, -gradk.component(1));
		Ak.replace(7,  gradk.component(0));
		Akhat = C_Ak*Ak;//Ak/(turbulence->epsilon()/sqrt(turbulence->k()));

		// omega gradient
		gradomega = fvc::grad(omega);

		// Other turbulence fields
		wallDistance = wallDist(mesh).y();

		// q's from ref Kaandorp 2020, 10.1016/j.compfluid.2020.104497
		q1 = 0.5 * (mag(R)*mag(R) - mag(S)*mag(S)) /(max(mag(S)*mag(S),SMALL_S2)); // Ratio of excess rotation rate to strain rate
		q2 = min((sqrt(k)*wallDistance/(50.0*turbulence->nu())),2.0);  // Wall distance based Reynolds Number
		q3 = T_t_ke * mag(S); // Ratio of turbulent time scale to mean strain time scale
		volSymmTensorField turbR = (2.0/3.0)*I*k - (nut)*dev(twoSymm(fvc::grad(U)));
		q4 = mag(turbR)/(k); //Ratio of total to 1/2 * normal Reynolds stresses (TKE)

		// Basis tensors
		B1 = Shat & Shat;
		B2 = Shat & Shat & Shat;
		B3 = Rhat & Rhat;
		B4 = Aphat & Aphat;
		B5 = Akhat & Akhat;
		B6 = Rhat & Rhat & Shat;
		B7 = Rhat & Rhat & Shat & Shat;
		B8 = Rhat & Rhat & Shat & Rhat & Shat & Shat;
		B9 = Aphat & Aphat & Shat;
		B10= Aphat & Aphat & Shat & Shat;
		B11= Aphat & Aphat & Shat & Aphat & Shat & Shat;
		B12= Akhat & Akhat & Shat;
		B13= Akhat & Akhat & Shat & Shat;
		B14= Akhat & Akhat & Shat & Akhat & Shat & Shat;
		B15= Rhat & Aphat;
		B16= Aphat & Akhat;
		B17= Rhat & Akhat;

		B18= Rhat & Aphat & Shat;
		B19= Rhat & Aphat & Shat & Shat;

		B20= Rhat & Rhat & Aphat & Shat;
		B21= Aphat & Aphat & Rhat & Shat;

		B22= Rhat & Rhat & Aphat & Shat & Shat;
		B23= Aphat & Aphat & Rhat & Shat & Shat;	

		B24= Rhat & Rhat & Shat & Aphat & Shat & Shat;
		B25= Aphat & Aphat & Shat & Rhat & Shat & Shat;

		B26= Rhat & Akhat & Shat;
		B27= Rhat & Akhat & Shat & Shat;

		B28= Rhat & Rhat & Akhat & Shat;
		B29= Akhat & Akhat & Rhat & Shat;

		B30= Rhat & Rhat & Akhat & Shat & Shat;
		B31= Akhat & Akhat & Rhat & Shat & Shat;	

		B32= Rhat & Rhat & Shat & Akhat & Shat & Shat;
		B33= Akhat & Akhat & Shat & Rhat & Shat & Shat;

		B34= Aphat & Akhat & Shat;
		B35= Aphat & Akhat & Shat & Shat;

		B36= Aphat & Aphat & Akhat & Shat;
		B37= Akhat & Akhat & Aphat & Shat;

		B38= Aphat & Aphat & Akhat & Shat & Shat;
		B39= Akhat & Akhat & Aphat & Shat & Shat;	

		B40= Aphat & Aphat & Shat & Akhat & Shat & Shat;
		B41= Akhat & Akhat & Shat & Aphat & Shat & Shat;

		B42= Rhat & Aphat & Akhat;
		
		B43= Rhat & Aphat & Akhat & Shat;
		B44= Rhat & Akhat & Aphat & Shat;
		B45= Rhat & Aphat & Akhat & Shat & Shat;
		B46= Rhat & Akhat & Aphat & Shat & Shat;
		B47= Rhat & Aphat & Shat & Akhat & Shat & Shat;

		// Invariants
		I1_1 = tr(B1);
		I1_2 = tr(B2);
		I1_3 = tr(B3);
		I1_4 = tr(B4);
		I1_5 = tr(B5);
		I1_6 = tr(B6);
		I1_7 = tr(B7);
		I1_8 = tr(B8);
		I1_9 = tr(B9);
		I1_10= tr(B10);
		I1_11= tr(B11);
		I1_12= tr(B12);
		I1_13= tr(B13);
		I1_14= tr(B14);
		I1_15= tr(B15);
		I1_16= tr(B16);
		I1_17= tr(B17);

		I1_18= tr(B18);
		I1_19= tr(B19);

		I1_20= tr(B20);
		I1_21= tr(B21);

		I1_22= tr(B22);
		I1_23= tr(B23);	

		I1_24= tr(B24);
		I1_25= tr(B25);

		I1_26= tr(B26);
		I1_27= tr(B27);

		I1_28= tr(B28);
		I1_29= tr(B29);

		I1_30= tr(B30);
		I1_31= tr(B31);	

		I1_32= tr(B32);
		I1_33= tr(B33);

		I1_34= tr(B34);
		I1_35= tr(B35);

		I1_36= tr(B36);
		I1_37= tr(B37);

		I1_38= tr(B38);
		I1_39= tr(B39);	

		I1_40= tr(B40);
		I1_41= tr(B41);

		I1_42= tr(B42);
		
		I1_43= tr(B43);
		I1_44= tr(B44);
		I1_45= tr(B45);
		I1_46= tr(B46);
		I1_47= tr(B47);

		

		I2_1 = 0.5*(tr(B1)*tr(B1) - tr (B1 & B1));
		I2_2 = 0.5*(tr(B2)*tr(B2) - tr (B2 & B2));
		I2_3 = 0.5*(tr(B3)*tr(B3) - tr (B3 & B3));
		I2_4 = 0.5*(tr(B4)*tr(B4) - tr (B4 & B4));
		I2_5 = 0.5*(tr(B5)*tr(B5) - tr (B5 & B5));
		I2_6 = 0.5*(tr(B6)*tr(B6) - tr (B6 & B6));
		I2_7 = 0.5*(tr(B7)*tr(B7) - tr (B7 & B7));
		I2_8 = 0.5*(tr(B8)*tr(B8) - tr (B8 & B8));
		I2_9 = 0.5*(tr(B9)*tr(B9) - tr (B9 & B9));
		I2_10 = 0.5*(tr(B10)*tr(B10) - tr (B10 & B10));
		I2_11 = 0.5*(tr(B11)*tr(B11) - tr (B11 & B11));
		I2_12 = 0.5*(tr(B12)*tr(B12) - tr (B12 & B12));
		I2_13 = 0.5*(tr(B13)*tr(B13) - tr (B13 & B13));
		I2_14 = 0.5*(tr(B14)*tr(B14) - tr (B14 & B14));
		I2_15 = 0.5*(tr(B15)*tr(B15) - tr (B15 & B15));
		I2_16 = 0.5*(tr(B16)*tr(B16) - tr (B16 & B16));
		I2_17 = 0.5*(tr(B17)*tr(B17) - tr (B17 & B17));
		I2_18 = 0.5*(tr(B18)*tr(B18) - tr (B18 & B18));
		I2_19 = 0.5*(tr(B19)*tr(B19) - tr (B19 & B19));
		I2_20 = 0.5*(tr(B20)*tr(B20) - tr (B20 & B20));
		I2_21 = 0.5*(tr(B21)*tr(B21) - tr (B21 & B21));
		I2_22 = 0.5*(tr(B22)*tr(B22) - tr (B22 & B22));
		I2_23 = 0.5*(tr(B23)*tr(B23) - tr (B23 & B23));
		I2_24 = 0.5*(tr(B24)*tr(B24) - tr (B24 & B24));
		I2_25 = 0.5*(tr(B25)*tr(B25) - tr (B25 & B25));
		I2_26 = 0.5*(tr(B26)*tr(B26) - tr (B26 & B26));
		I2_27 = 0.5*(tr(B27)*tr(B27) - tr (B27 & B27));
		I2_28 = 0.5*(tr(B28)*tr(B28) - tr (B28 & B28));
		I2_29 = 0.5*(tr(B29)*tr(B29) - tr (B29 & B29));
		I2_30 = 0.5*(tr(B30)*tr(B30) - tr (B30 & B30));
		I2_31 = 0.5*(tr(B31)*tr(B31) - tr (B31 & B31));
		I2_32 = 0.5*(tr(B32)*tr(B32) - tr (B32 & B32));
		I2_33 = 0.5*(tr(B33)*tr(B33) - tr (B33 & B33));
		I2_34 = 0.5*(tr(B34)*tr(B34) - tr (B34 & B34));
		I2_35 = 0.5*(tr(B35)*tr(B35) - tr (B35 & B35));
		I2_36 = 0.5*(tr(B36)*tr(B36) - tr (B36 & B36));
		I2_37 = 0.5*(tr(B37)*tr(B37) - tr (B37 & B37));
		I2_38 = 0.5*(tr(B38)*tr(B38) - tr (B38 & B38));
		I2_39 = 0.5*(tr(B39)*tr(B39) - tr (B39 & B39));
		I2_40 = 0.5*(tr(B40)*tr(B40) - tr (B40 & B40));
		I2_41 = 0.5*(tr(B41)*tr(B41) - tr (B41 & B41));
		I2_42 = 0.5*(tr(B42)*tr(B42) - tr (B42 & B42));
		I2_43 = 0.5*(tr(B43)*tr(B43) - tr (B43 & B43));
		I2_44 = 0.5*(tr(B44)*tr(B44) - tr (B44 & B44));
		I2_45 = 0.5*(tr(B45)*tr(B45) - tr (B45 & B45));
		I2_46 = 0.5*(tr(B46)*tr(B46) - tr (B46 & B46));
		I2_47 = 0.5*(tr(B47)*tr(B47) - tr (B47 & B47));
		/*I3_1 = det(B1);
		I3_2 = det(B2);
		I3_3 = det(B3);
		I3_4 = det(B4);
		I3_5 = det(B5);
		I3_6 = det(B6);
		I3_7 = det(B7);
		I3_8 = det(B8);
		I3_9 = det(B9);
		I3_10 = det(B10);
		I3_11 = det(B11);
		I3_12 = det(B12);
		I3_13 = det(B13);
		I3_14 = det(B14);
		I3_15 = det(B15);
		I3_16 = det(B16);
		I3_17 = det(B17);
		I3_18 = det(B18);
		I3_19 = det(B19);
		I3_20 = det(B20);
		I3_21 = det(B21);
		I3_22 = det(B22);
		I3_23 = det(B23);
		I3_24 = det(B24);
		I3_25 = det(B25);
		I3_26 = det(B26);
		I3_27 = det(B27);
		I3_28 = det(B28);
		I3_29 = det(B29);
		I3_30 = det(B30);
		I3_31 = det(B31);
		I3_32 = det(B32);
		I3_33 = det(B33);
		I3_34 = det(B34);
		I3_35 = det(B35);
		I3_36 = det(B36);
		I3_37 = det(B37);
		I3_38 = det(B38);
		I3_39 = det(B39);
		I3_40 = det(B40);
		I3_41 = det(B41);
		I3_42 = det(B42);
		I3_43 = det(B43);
		I3_44 = det(B44);
		I3_45 = det(B45);
		I3_46 = det(B46);
		I3_47 = det(B47);
		*/
		T1 = Shat;
		T2 = (Shat & Rhat) - (Rhat & Shat);
		T3 = (Shat & Shat) - onethird*(tr(Shat & Shat))*I;
		T4 = (Rhat & Rhat) - onethird*(tr(Rhat & Rhat))*I;
		T5 = (Rhat & Shat & Shat) - (Shat & Shat & Rhat);
		T6 = (Rhat & Rhat & Shat) + (Shat & Rhat & Rhat) - twothird*(tr(Shat & Rhat & Rhat))*I;
		T7 = (Rhat & Shat & Rhat & Rhat) - (Rhat & Rhat & Shat & Rhat);
		T8 = (Shat & Rhat & Shat & Shat) - (Shat & Shat & Rhat & Shat);
		T9 = (Rhat & Rhat & Shat & Shat) + (Shat & Shat & Rhat & Rhat) - twothird*(tr(Shat & Shat & Rhat & Rhat))*I;
		T10= (Rhat & Shat & Shat & Rhat & Rhat) - (Rhat & Rhat & Shat & Shat & Rhat);
		
		lambda1 = tr(Shat & Shat);
    	lambda2 = tr(Rhat & Rhat);
    	lambda3 = tr(Shat & Shat & Shat);
    	lambda4 = tr(Rhat & Rhat & Shat);
    	lambda5 = tr(Rhat & Rhat & Shat & Shat);

		Info<< "Writing fields....\n" << endl;
		T1.write();
		T2.write();
		T3.write();
		T4.write();
		T5.write();
		T6.write();
		T7.write();
		T8.write();
		T9.write();
		T10.write();
		/*B1.write();
		B2.write();
		B3.write();
		B4.write();
		B5.write();
		B6.write();
		B7.write();
		B8.write();
		B9.write();
		B10.write();
		B11.write();
		B12.write();
		B13.write();
		B14.write();
		B15.write();
		B16.write();
		B17.write();
		B18.write();
		B19.write();
		B20.write();
		B20.write();
		B21.write();
		B22.write();
		B23.write();
		B24.write();
		B25.write();
		B26.write();
		B27.write();
		B28.write();
		B29.write();
		B30.write();
		B31.write();
		B32.write();
		B33.write();
		B34.write();
		B35.write();
		B36.write();
		B37.write();
		B38.write();
		B39.write();
		B40.write();
		B41.write();
		B42.write();
		B43.write();
		B44.write();
		B45.write();
		B46.write();
		B47.write();*/
		I1_1.write();
		I1_2.write();
		I1_3.write();
		I1_4.write();
		I1_5.write();
		I1_6.write();
		I1_7.write();
		I1_8.write();
		I1_9.write();
		I1_10.write();
		I1_11.write();
		I1_12.write();
		I1_13.write();
		I1_14.write();
		I1_15.write();
		I1_16.write();
		I1_17.write();
		I1_18.write();
		I1_19.write();
		I1_20.write();
		I1_20.write();
		I1_21.write();
		I1_22.write();
		I1_23.write();
		I1_24.write();
		I1_25.write();
		I1_26.write();
		I1_27.write();
		I1_28.write();
		I1_29.write();
		I1_30.write();
		I1_31.write();
		I1_32.write();
		I1_33.write();
		I1_34.write();
		I1_35.write();
		I1_36.write();
		I1_37.write();
		I1_38.write();
		I1_39.write();
		I1_40.write();
		I1_41.write();
		I1_42.write();
		I1_43.write();
		I1_44.write();
		I1_45.write();
		I1_46.write();
		I1_47.write();
		I2_1.write();
		I2_2.write();
		I2_3.write();
		I2_4.write();
		I2_5.write();
		I2_6.write();
		I2_7.write();
		I2_8.write();
		I2_9.write();
		I2_10.write();
		I2_11.write();
		I2_12.write();
		I2_13.write();
		I2_14.write();
		I2_15.write();
		I2_16.write();
		I2_17.write();
		I2_18.write();
		I2_19.write();
		I2_20.write();
		I2_21.write();
		I2_22.write();
		I2_23.write();
		I2_24.write();
		I2_25.write();
		I2_26.write();
		I2_27.write();
		I2_28.write();
		I2_29.write();
		I2_30.write();
		I2_31.write();
		I2_32.write();
		I2_33.write();
		I2_34.write();
		I2_35.write();
		I2_36.write();
		I2_37.write();
		I2_38.write();
		I2_39.write();
		I2_40.write();
		I2_41.write();
		I2_42.write();
		I2_43.write();
		I2_44.write();
		I2_45.write();
		I2_46.write();
		I2_47.write();
		/*
		I3_1.write();
		I3_2.write();
		I3_3.write();
		I3_4.write();
		I3_5.write();
		I3_6.write();
		I3_7.write();
		I3_8.write();
		I3_9.write();
		I3_10.write();
		I3_11.write();
		I3_12.write();
		I3_13.write();
		I3_14.write();
		I3_15.write();
		I3_16.write();
		I3_17.write();
		I3_18.write();
		I3_19.write();
		I3_20.write();
		I3_21.write();
		I3_22.write();
		I3_23.write();
		I3_24.write();
		I3_25.write();
		I3_26.write();
		I3_27.write();
		I3_28.write();
		I3_29.write();
		I3_30.write();
		I3_31.write();
		I3_32.write();
		I3_33.write();
		I3_34.write();
		I3_35.write();
		I3_36.write();
		I3_37.write();
		I3_38.write();
		I3_39.write();
		I3_40.write();
		I3_41.write();
		I3_42.write();
		I3_43.write();
		I3_44.write();
		I3_45.write();
		I3_46.write();
		I3_47.write();
		*/

		q1.write();
		q2.write();
		q3.write();
		q4.write();
		S.write();
		Shat.write();
		Rhat.write();
		R.write();
		T_t_ke.write();
		T_t_nut.write();
		T_k.write();
		lambda1.write();
		lambda2.write();
		lambda3.write();
		lambda4.write();
		lambda5.write();
		gradp.write();
		gradp.write();
		gradk.write();
		gradU.write();
		gradomega.write();
		p.write();
		Ap.write();
		Ak.write();
		Aphat.write();
		Akhat.write();
		C_Ak.write();
		epsilon.write();
        runTime.write();
		DUDt.write();
		wallDistance.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
