    Info << "solving in subsurface domain ... " << endl;

    if(flowSwitch)
    {

     //read gradP from gradPDict, which is calculated from surface domain
     dimensionedVector gradP(gradPDict.lookup("gradP"));
    
     //debug
     Info << "gradP = " << gradP << endl;
   
     //multiplier for gradP. If gradP is switched on, then this variable = 1.0
     scalar gradPMultiplier=0.0;
     if(gradPSwitch) gradPMultiplier = 1.0;
    
     //predicted velocity with only gradP, not the pressure distribution yet.
     volVectorField HbyA
     (
         IOobject
         (
            "HbyA",
            runTime.timeName(),
            meshSubsurface,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         meshSubsurface,
         kSub & gradPMultiplier/muSub*rho_w*gradP 
     );

     //flow flux due to the predicted velocity
     surfaceScalarField phiHbyA
         (
            "phiHbyA",
            fvc::flux(HbyA)
         );

     //No need to constrain the pressure here as in interFoam. However, we should make sure
     //the gradP is in streamwise dirctioon, i.e. parallel to any wall.
     //Update the pressure BCs to ensure flux consistency
     //constrainPressure(p, U, phiHbyA, rAUf);

    // Non-orthogonal corrector loop
        while (simpleSubsurface.correctNonOrthogonal())
        {
         /* // fvScalarMatric pSubEqn       
          if (compressibilitySwitch)
            {
              fvScalarMatrix pSubEqn
              (
                   fvm::ddt(pSub)-fvm::laplacian(kSub/(muSub*rho_w)/(alpha-alpha*theta.value()+beta*theta.value()), pSub) == fvc::div(phiHbyA)
              );
              pSubEqn.solve();
            if (simpleSubsurface.finalNonOrthogonalIter())
           {
              //get the conservative flux
              //Here it should be just the flux returned by pSubEqn.flux().
              //Based on the pressure equation, the flux here is based on Darcy
              //velocity, not porewater velocity (which is what we want). So it should be 
              //divided by the porosity outside of the non-orthogonal corrector loop. 
              phiSub = phiHbyA + pSubEqn.flux();
           }

            }
 
           else
            {*/
             fvScalarMatrix pSubEqn
              (
                   fvm::laplacian(-kSub*k_Coeff/muSub*rho_w, pSub) == fvc::div(phiHbyA)//added the k_coeff for heterogeneity problem_BL it still works for 2021
              );
             pSubEqn.solve();   
            if (simpleSubsurface.finalNonOrthogonalIter())
           {
              //get the conservative flux
              //Here it should be just the flux returned by pSubEqn.flux().
              //Based on the pressure equation, the flux here is based on Darcy
              //velocity, not porewater velocity (which is what we want). So it should be 
              //divided by the porosity outside of the non-orthogonal corrector loop. 
              phiSub = phiHbyA + pSubEqn.flux();
           }

            }
/*       
           //pSubEqn.setReference(pRefCell, pRefValue);
           pSubEqn.solve();

           if (simpleSubsurface.finalNonOrthogonalIter())
           {
              //get the conservative flux
              //Here it should be just the flux returned by pSubEqn.flux().
              //Based on the pressure equation, the flux here is based on Darcy
              //velocity, not porewater velocity (which is what we want). So it should be 
              //divided by the porosity outside of the non-orthogonal corrector loop. 
              phiSub = phiHbyA + pSubEqn.flux();
           }
        }
*/
        //make the flux based on porewater velocity (needed in scalar transport equation).
        phiSub /= theta.value();
    
    
  

    //update porewater velocity (=Darcy velocity/porosity)
    //However, for incompressible surface flow solvers,
    //pressure is already divided by rho_water. So we 
    //need to multiply water density back
    if(gradPSwitch)
    {
     
    // gwform_heter_tensor USub == (-K & fvc::grad(pSub/g));
    // hyporheicFoam      USub == -kSub*k_Coeff/muSub*rho_w*(fvc::grad(pSub)-gradP)/theta.value();
//change for hetero
       //USub == -kSub & (fvc::grad(pSub)-gradP) * k_Coeff/muSub*rho_w/theta.value();
       USub == (-kSub & (fvc::grad(pSub)-gradP))* k_Coeff/muSub*rho_w/theta.value();
    }  
    else
    {  
//       USub == -kSub*k_Coeff/muSub*rho_w*fvc::grad(pSub)/theta.value();
       USub == (-kSub & fvc::grad(pSub)) * k_Coeff/muSub*rho_w/theta.value();

    }

 
//       USub == -kSub*k_Coeff/muSub*rho_w*(fvc::grad(pSub)-gradP)/theta.value();
  //     USub == -kSub*k_Coeff/muSub*rho_w*fvc::grad(pSub)/theta.value();
    

    Info<< "Continuity error of phiSub= "
        << mag(fvc::div(phiSub))().weightedAverage(meshSubsurface.V()).value()
        << endl;
}

    //Calculate Mechanical Dispersion Tensor referenced from Cardenas et al. 2011 Eq(4 and 5)
    tensor Kron=tensor::one;

    //a small velocity for division
    dimensionedScalar USmall
      (
         "USmall",
         dimensionSet(0,1,-1,0,0,0,0),
         VSMALL
      );

    //XL: not quite sure about this formula, regarding
    //    the definition of transverse and longitudinal.
    //    Eqn. (14) in Kaufman et al. (2017) WRR is different.
    //    Here, (USub * USub) is the outer product of velocity vector, whcih
    //    produces a symmetric tensor.
    volTensorField DSub_m("DSub_m",  Kron*(mag(USub)*alpha_T)
                                   + (alpha_L-alpha_T)/(mag(USub)+USmall)*(USub * USub));

    DSub_Coeff == DSub_m + Dm_tensor;

    //optional: kill the diffusive flux across the SWI
    if(ignoreDiffusiveFluxCouplingSwitch)
    {
        DSub_Coeff.boundaryFieldRef()[subSurfaceTopPatchID] *= 0.0;
    }

    // solve subSurface scalar transport equation
    if(CSubSwitch)
    { 
        #include "CSubEqn.H"

        //update scalar fluxes
        //this flux should be multiplied by the porosity because 
        //it is evaluated on the interface patch.
        advectiveCFluxSub = phiSub*fvc::interpolate(CSub)*theta.value();
    }
        // solve subSurface chemical reactive transport equation
    if(cheSubSwitch)
    {
        #include "YEqn.H"

        //update scalar fluxes
        //this flux should be multiplied by the porosity because 
        //it is evaluated on the interface patch.
   //     advectiveCFluxSub = phiSub*fvc::interpolate(CSub)*theta.value();
    }
