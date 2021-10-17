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

   //solve the unsteady, advection-diffiusion-reaction equationa// 20210112 remove the recycle of bags
    forAll(surY, i)
    {
        if (surcomposition.active(i))
        {
            volScalarField& surYi =surY[i];
            if (mag(runTime.value())>0)
            {
                fvScalarMatrix surYiEqn
                (

                    fvm::ddt(surYi)
              	  + fvm::div(phi, surYi)
                  - fvm::laplacian(Dt_Coeff, surYi)
                 ==
                    1*(surreaction->R(surYi))
                    + fvOptions(surYi)
		  
                );
     
                surYiEqn.relax();

            	fvOptions.constrain(surYiEqn);

            	surYiEqn.solve(meshSurface.solver("Yi"));
	
            	fvOptions.correct(surYi);

            	surYi.max(0.0);
            }
/*	    else if (0<=mag(runTime.value()) && mag(runTime.value())<0)
            
            {
		fvScalarMatrix surYiEqn
	        (
	             fvm::ddt(surYi)
	           + fvm::div(phi, surYi)
	           - fvm::laplacian(Dt_Coeff, surYi)
	           ==
	             TSur*(surreaction->R(surYi))
	           + fvOptions(surYi)
	         );
	   
           	 surYiEqn.relax();

           	 fvOptions.constrain(surYiEqn);
	
           	 surYiEqn.solve(meshSurface.solver("Yi"));

           	 fvOptions.correct(surYi);

           	 surYi.max(0.0);
	    }            
	    else if (0<=mag(runTime.value())&& mag(runTime.value())<0)
            {
                fvScalarMatrix surYiEqn
                (

                    fvm::ddt(surYi)
              	  + fvm::div(phi, surYi)
                  - fvm::laplacian(Dt_Coeff, surYi)
                 ==
                    TSur*(surreaction->R(surYi))
//                  + fvOptions(surYi)
		  
                );
     
                surYiEqn.relax();

  //          	fvOptions.constrain(surYiEqn);

            	surYiEqn.solve(meshSurface.solver("Yi"));
	
    //        	fvOptions.correct(surYi);

            	surYi.max(0.0);
            }
	    else if (0<=mag(runTime.value())&& mag(runTime.value())<0)
            
            {
		fvScalarMatrix surYiEqn
	        (
	             fvm::ddt(surYi)
	           + fvm::div(phi, surYi)
	           - fvm::laplacian(Dt_Coeff, surYi)
	           ==
	             surreaction->R(surYi)
	           + fvOptions(surYi)
	         );
	   
           	 surYiEqn.relax();

           	 fvOptions.constrain(surYiEqn);
	
           	 surYiEqn.solve(meshSurface.solver("Yi"));

           	 fvOptions.correct(surYi);

           	 surYi.max(0.0);
	    }   
	    else// if (92400<=mag(runTime.value())<870000)
            {
                fvScalarMatrix surYiEqn
                (

                    fvm::ddt(surYi)
              	  + fvm::div(phi, surYi)
                  - fvm::laplacian(Dt_Coeff, surYi)
                 ==
                    TSur*(surreaction->R(surYi)) // Tsur is the control of reactions 0/1
//                  + fvOptions(surYi)
		  
                );
     
                surYiEqn.relax();

  //          	fvOptions.constrain(surYiEqn);

            	surYiEqn.solve(meshSurface.solver("Yi"));
	
    //        	fvOptions.correct(surYi);

            	surYi.max(0.0);
            }
*/	    if(surYi.name()=="CO2")
	    {
	    	advectiveCO2FluxSurf = phi*fvc::interpolate(surYi);
	    }
            if(surYi.name()=="O2")
            {
                advectiveO2FluxSurf = phi*fvc::interpolate(surYi);
            }
			
        }
    }

