/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    hyporheicScalarFoam

Description
    A solver for hyporheic flow by combining pimpleFoam and groundWaterFoam
    1. pimpleFoam: Large time-step transient solver for incompressible, flow 
          using the PIMPLE (merged PISO-SIMPLE) algorithm. 
          Turbulence modelling, i.e. laminar, RAS or LES
    2. groundWaterFoam: Solves unsteady saturated porous media flow.

    Fluid is only one-way coupled, i.e., the subsurface flow is driven by 
    suface water pressure on SWI but the surface flow is not influenced by
    subsurface flow. In other words, the seepage at the bottom does not
    affect the surface flow.

    Scalar transport is added to the solver for both domains. The coupling 
    of scalar transport is through the boundary conditions:
    1. surfaceBottom
    2. subSurfaceTop

    Only one conservative scalar is solved. But this can be changed later 
    for cases like non-conservative (reactions) and multispecies. The porous
    medium is assumed to be isotropoic so the dispersivity tensor can be
    significantly simplified.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionModel.H"
#include "multivariateScheme.H"
#include "fvcSmooth.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#include "surfaceReaerationFvPatchScalarField.H"

#include "surfaceBottomFvPatchScalarField.H"
#include "subsurfaceTopFvPatchScalarField.H"
#include "chesurfaceBottomFvPatchScalarField.H"
#include "chesubsurfaceTopFvPatchScalarField.H"

#include "surfaceBottomFixedValueFvPatchScalarField.H"
#include "subsurfaceTopFixedValueFvPatchScalarField.H"
#include "chesurfaceBottomFixedValueFvPatchScalarField.H"
#include "chesubsurfaceTopFixedValueFvPatchScalarField.H"


#include "primitivePatchInterpolation.H"
#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"

    pimpleControl pimple(meshSurface);
    simpleControl simpleSubsurface(meshSubsurface);

    #include "createTimeControls.H"
    #include "createSurfaceFields.H"
    #include "createSubsurfaceFields.H"

    #include "createInterpolator.H"

    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    #include "readSimulationProperties.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //subiteration loop (used for coupling surface and subsurface domains)
        //fixed for now (need to add convergence check)
        for(int coupleIter=1;coupleIter<=maxCoupleIter;coupleIter++)
        {
            Info << nl <<  "Coupled iteration # " << coupleIter << " ..." << endl;

            //solve the surface flow equation
            #include "solve_surfaceFlow.H"

            //mapping bed bottom pressure and scalar from the surface 
            //flow region to the subsurface region
            #include "mappingSurfaceToBed.H"

            //solve the subsurface flow equation
            #include "solve_subSurfaceFlow.H"

            //mapping top scalar C from the subsurface flow region to 
            //the surface region
           #include "mappingBedToSurface.H"
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
