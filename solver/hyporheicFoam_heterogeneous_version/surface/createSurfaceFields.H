Info<< "Reading surface field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);
Info<< "Reading subsurface field p\n" << endl;

volScalarField TSur
(
    IOobject
    (
       "T",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);
Info<< "Reading surface field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

Info<< "Reading surface field C\n" << endl;
volScalarField CSurf //surface scalar 
(
    IOobject
    (
        "C",
        runTime.timeName(),
        meshSurface,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    meshSurface
);

#include "createSurfacePhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, meshSurface.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
meshSurface.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

//whether to solve scalar C transport in surface domain
Switch CSurfSwitch
(
    laminarTransport.lookup("CSurfSwitch")
);

//
//whether to solve chemical reaction transport in surface domain
Switch CheSurfSwitch
(
    laminarTransport.lookup("CheSurfSwitch")
);


Info << "Surface domain scalar transport switch: " << CSurfSwitch << endl;

// Laminar Prandtl number for scalar
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number for scalar
dimensionedScalar Prt("Prt", dimless, laminarTransport);

Info<< "Surface domain pressure gradient: " << endl;
//IOdictionary for storing/sharing mean pressure gradient between
//surface and subsurface
IOdictionary gradPDict
(
    IOobject
    (
        "gradPDict",
        runTime.constant(),
        meshSurface,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
);

dimensionedVector dummyGradP
(
     "dummyGradP",
     dimensionSet(0,1,-2,0,0,0,0),
     vector(0,0,0)
);

gradPDict.add("gradP",dummyGradP);


Info<< "Dt_Coeff...\n" << endl;
    volScalarField Dt_Coeff
    (
        IOobject
        (
            "Dt_Coeff",
            runTime.timeName(),
            meshSurface,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshSurface,
        dimensionedScalar("0", dimensionSet(0,2,-1,0,0,0,0),0.0)
    );


Info<< "Creating reaction model\n" << endl;

autoPtr<reactionModels::rhoReactionModel> surreaction
   (
        reactionModels::rhoReactionModel::New(meshSurface)
   );

        rhoReactionThermo& surthermo = surreaction->thermo();

        basicMultiComponentMixture& surcomposition = surthermo.composition();
        PtrList<volScalarField>& surY = surcomposition.Y();

// Set the turbulence into the reaction model
        surreaction->setTurbulence(turbulence());

        multivariateSurfaceInterpolationScheme<scalar>::fieldTable surfields;

        forAll(surY, i)
            {
               surfields.add(surY[i]);
            }
        PtrList<volScalarField> YiFromSub = surY;

        forAll(YiFromSub,i)
            {
               surfields.add(YiFromSub[i]);
               //surfields.add(test);
		// YiFromSub[i].rename(surY[i].name()+"FromSub");   
                   //}
               //YiFromSub[i] == surY[i].name()+"FromSub";    
            }				
