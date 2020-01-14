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

#include "Tuple2.H"

inline void Foam::InhibitionEfficiencies::checkCoeff() const
{
    forAll(*this, i)
    {
        if(operator[](i) < 0.0)
        {
            FatalErrorInFunction
                << "Negative inhibition constant error " << nl
                << "specie: " << species_[i] << nl
                << "inhibition constant:" << operator[](i) << nl
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::InhibitionEfficiencies::InhibitionEfficiencies
(
    const speciesTable& species,
    const scalarList& efficiencies
)
:
    scalarList(efficiencies),
    species_(species)
{
    if (size() != species_.size())
    {
        FatalErrorInFunction
            << "number of efficiencies = " << size()
            << " is not equal to the number of species " << species_.size()
            << exit(FatalError);
    }

    checkCoeff();
}


inline Foam::InhibitionEfficiencies::InhibitionEfficiencies
(
    const speciesTable& species,
    const dictionary& dict
)
:
    scalarList(species.size()),
    species_(species)
{
    if (dict.found("inhibition_coeffs"))
    {
        List<Tuple2<word, scalar>> coeffs(dict.lookup("inhibition_coeffs"));
        if (coeffs.size() != species_.size())
        {
            FatalErrorInFunction
                << "number of efficiencies = " << coeffs.size()
                << " is not equat to the number of species " << species_.size()
                << exit(FatalIOError);
        }

        forAll(coeffs, i)
        {
            operator[](species[coeffs[i].first()]) = coeffs[i].second();
        }
    }
    else
    {
        scalar defaultEff = readScalar(dict.lookup("defaultInhibitionCoeffs"));
        scalarList::operator=(defaultEff);
    }

    checkCoeff();
}


// * * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::InhibitionEfficiencies::M(const scalarList& c) const
{
    scalar M = 1.0;
    forAll(*this, i)
    {
        if(mag(operator[](i)) > VSMALL)
        {
           M *= (operator[](i)/(operator[](i)+c[i]));
        }
    }

    return M;
}


inline void Foam::InhibitionEfficiencies::write(Ostream& os) const
{
    List<Tuple2<word, scalar>> coeffs(species_.size());
    forAll(coeffs, i)
    {
        coeffs[i].first() = species_[i];
        coeffs[i].second() = operator[](i);
    }

    os.writeKeyword("inhibition_coeffs") << coeffs << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InhibitionEfficiencies& mnes
)
{
    mnes.write(os);
    return os;
}


// ************************************************************************* //
