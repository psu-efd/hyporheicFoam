/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::MonodInhibitionArrheniusReactionRate

Description
    Arrhenius reaction rate enhanced by Monod-type interation and inhibition.

SourceFiles
    MonodInhibitionArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef MonodInhibitionArrheniusReactionRate_H
#define MonodInhibitionArrheniusReactionRate_H

#include "ArrheniusReactionRate.H"
#include "MonodEfficiencies.H"
#include "InhibitionEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class MonodInhibitionArrheniusReactionRate;

Ostream& operator<<(Ostream&, const MonodInhibitionArrheniusReactionRate&);


/*---------------------------------------------------------------------------*\
              Class MonodInhibitionArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class MonodInhibitionArrheniusReactionRate
:
    ArrheniusReactionRate
{
    // Private data

        MonodEfficiencies MonodEfficiencies_;
        InhibitionEfficiencies InhibitionEfficiencies_;


public:

    // Constructors

        //- Construct from components
        inline MonodInhibitionArrheniusReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const MonodEfficiencies& mnes,
            const InhibitionEfficiencies& ines
        );

        //- Construct from dictionary
        inline MonodInhibitionArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "MonodInhibitionArrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const MonodInhibitionArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MonodInhibitionArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
