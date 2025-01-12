Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("",dimensionSet(0,2,-3,0,0,0,0),0.0)
);

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



volSymmTensorField S //Mean strain rate
(
	IOobject
	(
		"S",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("S",dimensionSet(0,0,-1,0,0,0,0),symmTensor::zero)
); 

volSymmTensorField turbR //Reynolds stress tensor
(
	IOobject
	(
		"turbR",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("turbR",dimensionSet(0,2,-2,0,0,0,0),symmTensor::zero)
); 


volTensorField R //Mean rotation rate
(
	IOobject
	(
		"R",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedTensor("R",dimensionSet(0,0,-1,0,0,0,0),tensor::zero)
); 

volTensorField gradU //Velocity gradient (note: transpose of the Jacobian is used in OpenFOAM)
(
	IOobject
	(
		"gradU",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedTensor("gradU",dimensionSet(0,0,-1,0,0,0,0),tensor::zero)
); 

volVectorField divturbR //Reynolds stress tensor divergence
(
	IOobject
	(
		"divturbR",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("divturbR",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
); 


volScalarField q1 // Used to be ratio of excess rotation to strain rate, but this field was very unstable.
(
    IOobject
    (
        "q1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("q1",dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField q2 // Wall distance based Reynolds number
(
    IOobject
    (
        "q2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("q2",dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField q3 // Ratio of turbulent time scale to mean strain time scale
(
    IOobject
    (
        "q3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("q3",dimensionSet(0,0,0,0,0,0,0),0.0)
);

volScalarField q4 // Ratio of total Reynolds stress to TKE
(
    IOobject
    (
        "q4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("q4",dimensionSet(0,0,0,0,0,0,0),0.0)
);





volScalarField T_t_ke
(
    IOobject
    (
        "T_t_ke",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
	dimensionedScalar("T_t_ke",dimensionSet(0,0,1,0,0,0,0),1.0)
);

volVectorField gradp 
(
	IOobject
	(
		"gradp",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("gradp",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
); 

volVectorField gradnut 
(
	IOobject
	(
		"gradnut",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("gradnut",dimensionSet(0,1,-1,0,0,0,0),vector::zero)
); 

volVectorField gradk 
(
	IOobject
	(
		"gradk",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("gradk",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
); 

volVectorField gradomega
(
	IOobject
	(
		"gradomega",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("gradomega",dimensionSet(0,-1,-1,0,0,0,0),vector::zero)
); 

volVectorField gradepsilon
(
	IOobject
	(
		"gradepsilon",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("gradepsilon",dimensionSet(0,1,-3,0,0,0,0),vector::zero)
); 



volScalarField wallDistance
(
	IOobject
	(
		"wallDistance",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedScalar("wallDistance",dimensionSet(0,1,0,0,0,0,0),0.0)
);


volVectorField DUDt
(
	IOobject
	(
		"DUDt",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedVector("DUDt",dimensionSet(0,1,-2,0,0,0,0),vector::zero)
);



dimensionedScalar SMALL_S2
(
	"small",
	dimensionSet(0,0,-2,0,0,0,0),
	scalar(1E-30)
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());



singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
