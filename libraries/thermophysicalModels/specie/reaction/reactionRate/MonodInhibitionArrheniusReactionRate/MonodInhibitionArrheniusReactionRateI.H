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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::MonodInhibitionArrheniusReactionRate::MonodInhibitionArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const MonodEfficiencies& mnes,
    const InhibitionEfficiencies& ines
)
:
    ArrheniusReactionRate(A, beta, Ta),
    MonodEfficiencies_(mnes),
    InhibitionEfficiencies_(ines)
{}


inline Foam::MonodInhibitionArrheniusReactionRate::MonodInhibitionArrheniusReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    ArrheniusReactionRate
    (
        species,
        dict
    ),
    MonodEfficiencies_(species, dict),
    InhibitionEfficiencies_(species, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::MonodInhibitionArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return
        MonodEfficiencies_.M(c)
       *InhibitionEfficiencies_.M(c)
       *ArrheniusReactionRate::operator()(p, T, c);
}


inline void Foam::MonodInhibitionArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    MonodEfficiencies_.write(os);
    InhibitionEfficiencies_.write(os);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MonodInhibitionArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
