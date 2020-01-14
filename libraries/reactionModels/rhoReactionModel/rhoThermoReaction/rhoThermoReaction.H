/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::rhoThermoReaction

Description
    Density-based thermo model wrapper for reaction models

SourceFiles
    rhoThermoReaction.C

\*---------------------------------------------------------------------------*/

#ifndef rhoThermoReaction_H
#define rhoThermoReaction_H

#include "autoPtr.H"
#include "rhoReactionModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                    class rhoThermoReaction Declaration
\*---------------------------------------------------------------------------*/

class rhoThermoReaction
:
    public rhoReactionModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoThermoReaction(const rhoThermoReaction&);

        //- Disallow default bitwise assignment
        void operator=(const rhoThermoReaction&);


protected:

    // Protected data

        //- Pointer to thermo model
        autoPtr<rhoReactionThermo> thermoPtr_;


public:

    // Constructors

        //- Construct from components
        rhoThermoReaction
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~rhoThermoReaction();


    // Member Functions

        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
