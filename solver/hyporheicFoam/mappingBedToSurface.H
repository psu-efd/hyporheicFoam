
        //map the subsurface top scalar to the surface bed bottom

        surfaceBottomFvPatchScalarField& tempCSurfp =
              dynamic_cast<surfaceBottomFvPatchScalarField&>
                      (CSurf.boundaryFieldRef()[surfaceBottomPatchID]);


        tempCSurfp.CFromSub() =
          bedToSurfaceInterpolate.faceInterpolate<scalar>
          (
             CSub.boundaryField()[subSurfaceTopPatchID]
          );
          
         //map the subsurface top nonconservatice scalar to the surface bed bottom
                forAll(surY, i)
          {
                if (surcomposition.active(i))
               {
                        //volScalarField& surYi =surY[i];
                   YiFromSub[i].boundaryFieldRef()[surfaceBottomPatchID] ==
                                bedToSurfaceInterpolate.faceInterpolate<scalar>
                                (
                                        Y[i].boundaryField()[subSurfaceTopPatchID]
                                );

                //volScalarField& Yi =Y[i];
                        chesurfaceBottomFvPatchScalarField& tempYiSurfp =
                    dynamic_cast<chesurfaceBottomFvPatchScalarField&>(surY[i].boundaryFieldRef()[surfaceBottomPatchID]);

                        tempYiSurfp.YiFromSub() =
                                bedToSurfaceInterpolate.faceInterpolate<scalar>
                                (
                                        Y[i].boundaryField()[subSurfaceTopPatchID]
                                );
              }
          }


        const surfaceVectorField Sn_sub(meshSubsurface.Sf()/meshSubsurface.magSf());
        const vectorField& Sn_sub_f = Sn_sub.boundaryField()[subSurfaceTopPatchID];

        //map the subsurface top velocity to the surface bed bottom
        //Using Darcy velocity for consistency in fluid flux across 
        //the interface. Since USub is the pore-water velocity,
        //it needs to be multiplied by porosity to get Darcy velocity.  
        //May not work if the upwelling/downwelling is significant!
        //Only use the normla component may help because the only thing matters
        //is the normal velocity for flux.

      if(velocityNormalCouplingSwitch)
        {
            //Use the patch normal component.
            U.boundaryFieldRef()[surfaceBottomPatchID] ==
               bedToSurfaceInterpolate.faceInterpolate<vector>
               (
                  (Sn_sub_f & USub.boundaryField()[subSurfaceTopPatchID]*theta.value())*Sn_sub_f
               );
        }
        else
        {
            //Use the whole vector.
            U.boundaryFieldRef()[surfaceBottomPatchID] ==
               bedToSurfaceInterpolate.faceInterpolate<vector>
               (
                  USub.boundaryField()[subSurfaceTopPatchID]*theta.value()
               );
        }
       
        //Info << "U.boundaryField = " << U.boundaryField()[surfaceBottomPatchID] << endl; 
