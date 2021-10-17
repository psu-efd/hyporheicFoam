    Info<< "Reading simulationProperties\n" << endl;

    dictionary dict1(IFstream("constant/simulationProperties")());

    //whether to solve flow
    //on: flow will be simulated.
    //off: flow will not be simulated (an option for using pre-simulated flow fields)
    Switch flowSwitch
    (
        dict1.lookupOrDefault("flowSwitch",true)
    );

    Info << "flowSwitch = " << flowSwitch << endl;

    if(flowSwitch)
    {
        Info << "Flow (in both surface and subsurface domains) will be simulated." << endl;
    }
    else
    {
        Info << "Flow (in both surface and subsurface domains) will not be simulated." << endl;
    }
 
    Switch compressibilitySwitch
    (
         dict1.lookupOrDefault("compressibilitySwitch",true)
    );
    Info << "compressibilitySwitch = " << compressibilitySwitch << endl;
  
    if(compressibilitySwitch)
       {
          Info << "compressibility is considered. " << endl;
       }
     else
       {
        Info << "compressibility is not considered. " << endl;
       }


 
    //whether to solve chemical reaction
    //On: chemical reaction will be simulated.
    //Off: chemical reaction will not be simulated.
    Switch cheSwitch
    (
        dict1.lookupOrDefault("cheSwitch",true)
    );
     
     Info << "cheSwitch = " << cheSwitch << endl;
     
     if(cheSwitch)
     {
         Info << "Chemical (in both surface and subsurface domains) will be simulated." << endl;
     }
     else
     {
         Info << "Chemical (in both surface and subsurface domains) will not be simulated." << endl;
     }
 


    //maxCoupleIter: maximum subiteration number
    //The default is 1 (just one iteration)
    label maxCoupleIter
    (
       dict1.lookupOrDefault("maxCoupleIter",1)
    );

    Info << "maxCoupleIter = " << maxCoupleIter << endl;

    //For the flow coupling, whether to only use the patch normal
    //component to the velocity. The default is true.
    Switch velocityNormalCouplingSwitch
    (
        dict1.lookupOrDefault("velocityNormalCouplingSwitch",true)
    );

    Info << "velocityNormalCouplingSwitch = " << velocityNormalCouplingSwitch << endl;

    if(velocityNormalCouplingSwitch)
    {
        Info << "Patch normal velocity is used for coupling at interface." << endl;
    }
    else
    {
        Info << "Full velocity vector is used for coupling at interface." << endl;
    }



    //For the scalar coupling, whether to kill the diffusive flux
    //The default is true.
    Switch ignoreDiffusiveFluxCouplingSwitch
    (
        dict1.lookupOrDefault("ignoreDiffusiveFluxCouplingSwitch",true)
    );

    Info << "ignoreDiffusiveFluxCouplingSwitch = " << ignoreDiffusiveFluxCouplingSwitch << endl;

    if(ignoreDiffusiveFluxCouplingSwitch)
    {
        Info << "Diffusive flux is ignored at interface." << endl;
    }
    else
    {
        Info << "Diffusive flux is not ignored at interface." << endl;
    }

