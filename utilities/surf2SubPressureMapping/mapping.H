
        Info<< "Mapping pressure fields from source time "
            << runTimeSource.timeName() << endl;

        Info<< "Reading surface field p\n" << endl;
        IOobject pSurfHeader
        (
            "p",
            runTimeSource.timeName(),
            meshSource,
            IOobject::MUST_READ
        );

        if (!pSurfHeader.typeHeaderOk<volScalarField>(true))
        {
            Info<< "    No surface pressure field" << endl;
            FatalError.exit();
        }

        volScalarField pSurf
        (
            pSurfHeader,
            meshSource
        );

        Info<< "Reading subsurface field p\n" << endl;
        volScalarField pSub
        (
            IOobject
            (
               "p",
               "0",   //always read subsurface at time = 0
               meshTarget,
               IOobject::MUST_READ,
               IOobject::AUTO_WRITE
            ),
            meshTarget
       );

       //mapping bed bottom pressure from the surface 
       //flow region to the subsurface region
       #include "mappingSurfaceToBed.H"

       pSub.write();
