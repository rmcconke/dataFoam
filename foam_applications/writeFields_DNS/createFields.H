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

// Reynolds stress tensor
Info<< "Reading field TauDNS\n" << endl;
volSymmTensorField tau 
(
	IOobject
	(
		"TauDNS",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
),
	mesh
); 

volSymmTensorField a
(
	IOobject
	(
		"a",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("a",dimensionSet(0,2,-2,0,0,0,0),symmTensor::zero)
); 

volSymmTensorField b
(
	IOobject
	(
		"b",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("b",dimensionSet(0,0,0,0,0,0,0),symmTensor::zero)
); 

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("k",dimensionSet(0,2,-2,0,0,0,0),0.0)
);

volTensorField gradU //Mean velocity gradient
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
