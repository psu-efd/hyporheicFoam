    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            meshSubsurface,
            fields,
            phiSub,
            meshSubsurface.divScheme("div(phiSub,Yi)")
        )
   );
 //update the reaction rate
    reaction->correct();


//label patchID1 = meshSubsurface.boundaryMesh().findPatchID("left");
//label patchID2 = meshSubsurface.boundaryMesh().findPatchID("right");
//Info << " patchID1111111 " << patchID1 <<endl;
//Info << " patchID2222222 " << patchID2 <<endl;
while (simpleSubsurface.correctNonOrthogonal())
{
    forAll(Y, i)
    {
        if (composition.active(i))
        {
            volScalarField& Yi = Y[i];
            if(Yi.name()=="CH2OO")
            {
                Yi.primitiveFieldRef() = 1;
                //forAll(Yi.boundaryFieldRef(),patchi)
		forAll(Yi.boundaryFieldRef(),patchi)                
		{
                    Yi.boundaryFieldRef()[patchi] = 100000;
                }
                //Yi.primitiveFieldRef() = 200;
            }
            else if (Yi.name()=="CH2O")
            {
                fvScalarMatrix YiEqn
                (

                   fvm::ddt(Yi)
                 + fvm::div(phiSub, Yi)
                 - fvm::laplacian(DSub_Coeff, Yi)
                 ==
                   reaction->R(Yi)//from the reaction
                 //  + fvOptions(Yi)//increasing in the specific zone.
		  //+ sourceValue // increasing in the whole domain.
                );
     
                YiEqn.relax();

               // fvOptions.constrain(YiEqn);

                YiEqn.solve(meshSubsurface.solver("Yi"));

//                fvOptions.correct(Yi);

                Yi.max(0.0);
            }else
            {
                fvScalarMatrix YiEqn
                (

                  fvm::ddt(Yi)
                 + fvm::div(phiSub, Yi)
                  - fvm::laplacian(DSub_Coeff, Yi)
                 ==
                    reaction->R(Yi)
                    //+ fvOptions(Yi)
                );
     
                YiEqn.relax();

                //fvOptions.constrain(YiEqn);

                YiEqn.solve(meshSubsurface.solver("Yi"));

               // fvOptions.correct(Yi);

                Yi.max(0.0);
            }
	    if(Yi.name()=="CO2")
            {
                advectiveCO2FluxSub = phiSub*fvc::interpolate(Yi);
            }
            if(Yi.name()=="O2")
            {
                advectiveO2FluxSub = phiSub*fvc::interpolate(Yi);
            }
        }
     }  
} 

