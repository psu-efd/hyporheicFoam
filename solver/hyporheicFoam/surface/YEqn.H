tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            meshSurface,
            fields,
            phi,
            meshSurface.divScheme("div(phi,Yi)")
        )
    );
    surreaction->correct();

   //solve the unsteady, advection-diffiusion-reaction equation
    forAll(surY, i)
    {
        if (surcomposition.active(i))
        {
            volScalarField& surYi =surY[i];
            
            fvScalarMatrix YiEqn
            (
                fvm::ddt(surYi)
              + fvm::div(phi, surYi)
              - fvm::laplacian(Dt_Coeff, surYi)
             ==
                surreaction->R(surYi)
              + fvOptions(surYi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(meshSurface.solver("Yi"));

            fvOptions.correct(surYi);

            surYi.max(0.0);
	    if(surYi.name()=="CO2")
	    {
	    	advectiveCO2FluxSurf = phi*fvc::interpolate(surYi);
	    }
            if(surYi.name()=="O2")
            {
                advectiveO2FluxSurf = phi*fvc::interpolate(surYi);
            }
			
        }
    }

