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




while (simpleSubsurface.correctNonOrthogonal())
{
    forAll(Y, i)
    {
        if (composition.active(i))
        {
            volScalarField& Yi = Y[i];
            if(Yi.name()=="CH2O")
            {
                Yi.primitiveFieldRef() = 1;
                forAll(Yi.boundaryFieldRef(),patchi)
                {
                    Yi.boundaryFieldRef()[patchi] = 1;
                }
                //Yi.primitiveFieldRef() = 200;
            }
            else
            {
                fvScalarMatrix YiEqn
                (

                  fvm::ddt(Yi)
                 + fvm::div(phiSub, Yi)
                  - fvm::laplacian(DSub_Coeff, Yi)
                 ==
                    reaction->R(Yi)
                  + fvOptions(Yi)
                );
     
                YiEqn.relax();

                fvOptions.constrain(YiEqn);

                YiEqn.solve(meshSubsurface.solver("Yi"));

                fvOptions.correct(Yi);

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

