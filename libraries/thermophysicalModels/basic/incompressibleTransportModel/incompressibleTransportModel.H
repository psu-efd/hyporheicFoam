/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

Class
    Foam::incompressibleTransportModel

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    incompressibleTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTransportModel_H
#define incompressibleTransportModel_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleTransportModel Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTransportModel
{
    // Private Member Functions

        //- Disallow copy construct
        incompressibleTransportModel(const incompressibleTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const incompressibleTransportModel&);


public:

    //- Runtime type information
    TypeName("incompressibleTransportModel");


    // Constructors

        //- Construct from components
        incompressibleTransportModel
        ();


    //- Destructor
    virtual ~incompressibleTransportModel();


    // Member Functions

        //- Return the dynamic laminar viscosity
        virtual tmp<volScalarField> mu() const = 0;

        //- Return the dynamic laminar viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const = 0;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
