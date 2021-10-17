    Info<< "Reading subsurface field p\n" << endl;

    volScalarField pSub
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            meshSubsurface,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshSubsurface
    );

    //porewater velocity (=Darcy velocity/porosity)
    volVectorField USub
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            meshSubsurface,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshSubsurface
    );

    volScalarField CSub  //Scalar
    (
        IOobject
        (
            "C",
            runTime.timeName(),
            meshSubsurface,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshSubsurface
    );
    volScalarField k_Coeff  //Scalar
    (
        IOobject
        (
            "k_Coeff",
            runTime.timeName(),
            meshSubsurface,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        meshSubsurface
    );
    #include "createSubsurfacePhi.H"

    Info<< "Reading subsurface transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            meshSubsurface,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    //whether to solve scalar C transport in subsurface domain
    Switch CSubSwitch 
    (
        transportProperties.lookup("CSubSwitch")
    );

     //whether to solve chemical reaction transport in subsurface domain
    Switch cheSubSwitch
    (
        transportProperties.lookup("cheSubSwitch")
    );


    //whether to apply the surface pressure gradient in subsurface domain
    Switch gradPSwitch 
    (
        transportProperties.lookup("gradPSwitch")
    );




   Info<< "Creating reaction model\n" << endl;

      autoPtr<reactionModels::rhoReactionModel> reaction
      (
         reactionModels::rhoReactionModel::New(meshSubsurface)
      );

      rhoReactionThermo& thermo = reaction->thermo();

      basicMultiComponentMixture& composition = thermo.composition();
      PtrList<volScalarField>& Y = composition.Y();



    Info << "Subsurface domain scalar transport switch: " << CSubSwitch << endl;


     
    //intrinsic permeability (L^2)
    dimensionedTensor kSub //change the dimensionedScalar K to dimensionedTensor K for the heterogeneous problems. 20210502 _Bing
    (
        transportProperties.lookup("k")
    );
    dimensionedScalar sourceValue
    (
        transportProperties.lookup("sourceValue")
    );
    //sediment matrix compressibility (1/Pressure)
    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );

    // porewater compressibility (1/Pressure)
    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    );

    // porosity
    dimensionedScalar theta
    (
        transportProperties.lookup("theta")
    );

    // porewater dynamic viscosity (Pressure*time)
    dimensionedScalar muSub
    (
        transportProperties.lookup("mu")
    );

    // define some constants
    // water density
    dimensionedScalar rho_w
    (
       "rho_w",
       dimensionSet(1,-3,0,0,0,0,0),
       scalar(1000.0)
    );

    // gravity constant
    dimensionedScalar g
    (
        "g",
        dimensionSet(0,1,-2,0,0,0,0),
        scalar(9.81)
    );

    // Calculate the hydraulic conductivity
    //dimensionedScalar K(kSub*rho_w*g/muSub);
   // Info << "Calculated hydraulic conductivity K = " << K << endl;

    // longitudinal dispersivity (Length)
    dimensionedScalar alpha_L
    (
        transportProperties.lookup("alpha_L")
    );

    // transverse dispersivity (Length)
    dimensionedScalar alpha_T
    (
        transportProperties.lookup("alpha_T")
    );

    // pore-water molecular diffusion coefficient
    dimensionedScalar Dm
    (
        transportProperties.lookup("Dm")
    );
  
    //molecular diffusion tensor
    dimensionedTensor Dm_tensor
        (
            "Dm_tensor",
            dimensionSet(0,2,-1,0,0,0,0),
            tensor(I)*Dm.value()
        );

    volTensorField DSub_Coeff 
    (
        IOobject
        (
            "DSub_Coeff",
            runTime.timeName(),
            meshSubsurface,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        meshSubsurface,
        dimensionedTensor("0", dimensionSet(0,2,-1,0,0,0,0),tensor::zero) 
    );


     //Info<< "Creating turbulence model\n" << endl;
     //autoPtr<incompressible::turbulenceModel> turbulence
     //(
     //  incompressible::turbulenceModel::New(USub, phiSub, laminarTransport)
     //);

     // Set the turbulence into the reaction model
      reaction->setTurbulence(turbulence());

      multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

      forAll(Y, i)
       {
         fields.add(Y[i]);
       }
      PtrList<volScalarField> YiFromSurf = Y;
      forAll(YiFromSurf, i)
       {
      fields.add(YiFromSurf[i]);
       }

