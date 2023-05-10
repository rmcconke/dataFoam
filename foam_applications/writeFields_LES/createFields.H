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

Info<< "Reading field UMean\n" << endl;
volVectorField UMean
(
    IOobject
    (
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field kMean_model\n" << endl;
volScalarField kMean_model
(
    IOobject
    (
        "kMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField kMean_tauMean
(
    IOobject
    (
        "kMean_tauMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kMean_tauMean",dimensionSet(0,2,-2,0,0,0,0),0.0)
);


Info<< "Reading field nutMean\n" << endl;
volScalarField nutMean
(
    IOobject
    (
        "nut", //Note: needs to be changed once we have updated LES data
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volSymmTensorField tauMean
(
	IOobject
	(
		"tauMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("tauMean",dimensionSet(0,2,-2,0,0,0,0),symmTensor::zero)
); 

volSymmTensorField aMean
(
	IOobject
	(
		"aMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("aMean",dimensionSet(0,2,-2,0,0,0,0),symmTensor::zero)
); 

volSymmTensorField subgrid_tauMean
(
	IOobject
	(
		"subgrid_tauMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("subgrid_tauMean",dimensionSet(0,2,-2,0,0,0,0),symmTensor::zero)
); 

volSymmTensorField bMean
(
	IOobject
	(
		"bMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("bMean",dimensionSet(0,0,0,0,0,0,0),symmTensor::zero)
); 

volTensorField gradUMean //Mean velocity gradient
(
	IOobject
	(
		"gradUMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedTensor("gradUMean",dimensionSet(0,0,-1,0,0,0,0),tensor::zero)
); 

volTensorField RMean //Mean rotation rate
(
	IOobject
	(
		"RMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedTensor("RMean",dimensionSet(0,0,-1,0,0,0,0),tensor::zero)
); 

volSymmTensorField SMean //Mean strain rate
(
	IOobject
	(
		"SMean",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
),
	mesh,
	dimensionedSymmTensor("SMean",dimensionSet(0,0,-1,0,0,0,0),symmTensor::zero)
); 

Info<< "Reading field UPrime2Mean\n" << endl;
volSymmTensorField UPrime2Mean
(
    IOobject
    (
        "UPrime2Mean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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
