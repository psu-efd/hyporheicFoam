/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    createPhi

Description
    Creates and initialises the relative face-flux field phi.

\*---------------------------------------------------------------------------*/

#ifndef createPhi_H
#define createPhi_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Info<< "Reading/calculating surface face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        meshSurface,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & meshSurface.Sf()
);
Info<< "Reading/calculating surface face flux field advectiveCO2FluxSurf\n" << endl;

surfaceScalarField advectiveCO2FluxSurf
(
        IOobject
        (
                "advectiveCO2FluxSurf",
                runTime.timeName(),
                meshSurface,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
       (linearInterpolate(U) & meshSurface.Sf())
);
Info<< "Reading/calculating surface face flux field advectiveO2FluxSurf\n" << endl;

surfaceScalarField advectiveO2FluxSurf
(
        IOobject
        (
                "advectiveO2FluxSurf",
                runTime.timeName(),
                meshSurface,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
        ),
       (linearInterpolate(U) & meshSurface.Sf())
);


//scalar advective flux (=C*phi) for the purpose of postprocessing
surfaceScalarField advectiveCFluxSurf
(
    IOobject
    (
        "advectiveCFluxSurf",
        runTime.timeName(),
        meshSurface,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (linearInterpolate(U) & meshSurface.Sf())*linearInterpolate(CSurf)
);

//scalar diffusive flux (=D*grad(C)*Sf) for the purpose of postprocessing
surfaceScalarField diffusiveCFluxSurf
(
    IOobject
    (
        "diffusiveCFluxSurf",
        runTime.timeName(),
        meshSurface,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*(linearInterpolate(U) & meshSurface.Sf())*linearInterpolate(CSurf)
);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
