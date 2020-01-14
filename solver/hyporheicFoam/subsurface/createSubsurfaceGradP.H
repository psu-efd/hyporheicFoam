    Info<< "Reading surface gradientProperties\n" << endl;
    IOdictionary gradientProperties
    (
        IOobject
        (
            "gradientProperties",
            runTime.constant(),
            meshSurface,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //read gradP from constant/surface/gradientProperties
    dimensionedVector gradP_
    (
       gradientProperties.lookup("gradP")
    );
