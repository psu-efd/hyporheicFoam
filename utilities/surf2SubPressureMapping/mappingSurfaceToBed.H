

  pSub.boundaryFieldRef()[subSurfaceTopPatchID] ==
             surfaceToBedInterpolate.faceInterpolate<scalar>
             (
                pSurf.boundaryField()[surfaceBottomPatchID] 
             );


