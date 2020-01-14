        Info << "solving in surface domain ... " << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if(flowSwitch)
            { 
                #include "UEqn.H"

               // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }

                //effective diffusion coefficient for scalar
                Dt_Coeff == turbulence->nu()/Pr + turbulence->nut()/Prt;

         
                //optional: kill the diffusive flux across the SWI
                if(ignoreDiffusiveFluxCouplingSwitch)
                {
                    Dt_Coeff.boundaryFieldRef()[surfaceBottomPatchID] *= 0.0;
                }
            }

            // scalar transport equation
            if(CSurfSwitch)
            {
                #include "CEqn.H"

                //update scalar fluxes
                advectiveCFluxSurf = phi*fvc::interpolate(CSurf);
            }



            if(CheSurfSwitch)
            {
                #include "YEqn.H"

                //update scalar fluxes
               // test = phi*fvc::interpolate(CO2);
            }

        }
