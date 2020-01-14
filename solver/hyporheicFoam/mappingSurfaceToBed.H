
        //map the surface bottom pressure to the sediment region top patch
        if(flowSwitch) 
        {
           pSub.boundaryFieldRef()[subSurfaceTopPatchID] ==
             surfaceToBedInterpolate.faceInterpolate<scalar>
             (
                p.boundaryField()[surfaceBottomPatchID] 
             );
        }


        //map the surface bottom scalar to the sediment region top patch
        subsurfaceTopFvPatchScalarField& tempCSubp =
              dynamic_cast<subsurfaceTopFvPatchScalarField&>(CSub.boundaryFieldRef()[subSurfaceTopPatchID]);

        tempCSubp.CFromSurf() =
            surfaceToBedInterpolate.faceInterpolate<scalar>
              (
                 CSurf.boundaryField()[surfaceBottomPatchID]
              );
   
         //map the surface bottom nonconservative scalar to the sediment region top patch
        {
           forAll(Y, i)
             {
                if (composition.active(i))
                {
                   //volScalarField& subYi =Y[i];

                   YiFromSurf[i].boundaryFieldRef()[subSurfaceTopPatchID] ==
                         surfaceToBedInterpolate.faceInterpolate<scalar>
                         (
                         surY[i].boundaryField()[surfaceBottomPatchID]
                         );
                   //forAll(surY, i)
                 //{ 
               //if(surcomposition.active(i))
            //{       
              //volScalarField& surYi =surY[i]; 
              chesubsurfaceTopFvPatchScalarField& tempsubYip =
                  dynamic_cast<chesubsurfaceTopFvPatchScalarField&>(Y[i].boundaryFieldRef()[subSurfaceTopPatchID]);

              tempsubYip.YiFromSurf() =
                 surfaceToBedInterpolate.faceInterpolate<scalar>
                 (
                  surY[i].boundaryField()[surfaceBottomPatchID]
                 );
               }
             }
        }

//   Info << "CSurf.boundaryField()[surfaceBottomPatchID] = " << CSurf.boundaryField()[surfaceBottomPatchID] << endl;
//   Info << "tempCSubp.CFromSurf = " << tempCSubp.CFromSurf() << endl;    
//   Info << "interpolate return = " <<  surfaceToBedInterpolate.faceInterpolate<scalar>(CSurf.boundaryField()[surfaceBottomPatchID]) << endl;



