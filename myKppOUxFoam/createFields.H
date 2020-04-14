Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    // dimension(0,0,0,0,0,0,0),
    mesh
);

volScalarField xi
(
    IOobject
    (
        "xi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField velInit
(
    IOobject
    (
        "velInit",
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

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookup("DT")
);

Info<< "Reading diffusivity Dxi \n" << endl;

dimensionedScalar Dxi
(
    transportProperties.lookup("Dxi")
);

Info<< "Reading chemical rate constant DK\n" << endl;

dimensionedScalar DK
(
    transportProperties.lookup("DK")
);

Info<< "Reading chemical rate constant Db\n" << endl;

dimensionedScalar Db
(
    transportProperties.lookup("Db")
);

Info<< "Reading rate repr dtheta\n" << endl;

dimensionedScalar dtheta
(
    transportProperties.lookup("dtheta")
);

Info<< "Reading mean stochastic value dmean\n" << endl;

dimensionedScalar dmean
(
    transportProperties.lookupOrDefault("dmean",0.0)
);

Info<< "Reading mean stochastic value dsigma\n" << endl;

dimensionedScalar dsigma
(
    transportProperties.lookup("dsigma")
);

Info<< "Reading mean velocity value barVel \n" << endl;

dimensionedScalar barVel
(
    transportProperties.lookup("barVel")
);
Info<< "Reading mean velocity value barVel \n" << endl;

dimensionedScalar vShift
(
    transportProperties.lookup("vShift")
);

// Set the option to choose the model of Stochastic field generation
//
word sdeDef("Ornstein");
word sdeScheme
(
    transportProperties.lookupOrDefault("sdeScheme",sdeDef )
);

dimensionedScalar dtau
(
    transportProperties.lookupOrDefault("dtau", dtheta)
);


bool useOldField =
    runTime.controlDict().lookupOrDefault("useOldField", false);


#include "createPhi.H"

#include "createFvOptions.H"


volScalarField volDivPhi
(
    IOobject
    (
        "volDivPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(phi)
    // mesh,
    // dimensionedScalar("volDivPhi", dimless, scalar(0.0))
);

volScalarField volDivU
(
    IOobject
    (
        "volDivU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(U)
    // mesh,
    // dimensionedScalar("volDivU", dimless, scalar(0.0))
);
