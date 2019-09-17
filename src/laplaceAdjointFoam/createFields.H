Info << "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
	"physicalProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

dimensionedScalar k(physicalProperties.lookup("k"));
dimensionedScalar beta(physicalProperties.lookup("beta"));
dimensionedScalar lambda(physicalProperties.lookup("lambda"));


/***************************************************************************/

Info << "Reading numericalParameters\n" << endl;
IOdictionary numericalParameters
(
    IOobject
    (
	"numericalParameters",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

scalar tol(numericalParameters.lookupOrDefault<scalar>("tol", 1e-8));
scalar alpha(numericalParameters.lookupOrDefault<scalar>("alpha", 1));
scalar c1(numericalParameters.lookupOrDefault<scalar>("c1", 1e-4));
scalar c2(numericalParameters.lookupOrDefault<scalar>("c2", 0.8));
/****************************************************************************/

Info << "Reading target field\n" << endl;
volScalarField yd
(
    IOobject
    (
	"yd",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading source field\n" << endl;
volScalarField f
(
    IOobject
    (
	"f",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);

Info << "Reading uMax field\n" << endl;
volScalarField uMax
(
    IOobject
    (
	"uMax",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);


Info << "Reading uMin field\n" << endl;
volScalarField uMin
(
    IOobject
    (
	"uMin",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);


/****************************************************************************/

Info << "Reading field y\n" << endl;
volScalarField y
(
    IOobject
    (
	"y",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading adjoint field\n" << endl;
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

Info << "Reading control u\n" << endl;
volScalarField u
(
    IOobject
    (
	"u",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

scalarField volField = mesh.V();

volScalarField ud
(
    IOobject
    (
	"ud",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField pd
(
    IOobject
    (
	"pd",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField uDiff
(
    IOobject
    (
	"uDiff",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   mag( u - ud)
);

volScalarField pDiff
(
    IOobject
    (
	"pDiff",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   mag( p - pd)
);

volScalarField yDiff
(
    IOobject
    (
	"yDiff",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
   mag( y - yd)
);


volScalarField uk
(
    IOobject
    (
	"uk",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    u
);

wordList patchesToIntegrate(physicalProperties.lookup("boundaryPatchList"));

