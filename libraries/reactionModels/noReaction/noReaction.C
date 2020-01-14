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

#include "noReaction.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactThermoType>
Foam::reactionModels::noReaction<ReactThermoType>::noReaction
(
    const word& modelType,
    const fvMesh& mesh,
    const word& reactionProperties,
    const word& phaseName
)
:
    ReactThermoType(modelType, mesh, phaseName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactThermoType>
Foam::reactionModels::noReaction<ReactThermoType>::~noReaction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactThermoType>
void Foam::reactionModels::noReaction<ReactThermoType>::correct()
{}


template<class ReactThermoType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionModels::noReaction<ReactThermoType>::R
(
    volScalarField& Y
) const
{
    tmp<fvScalarMatrix> tSu
    (
        new fvScalarMatrix(Y, dimless/dimTime)
        //new fvScalarMatrix(Y, dimMass/dimTime)
    );

    return tSu;
}


template<class ReactThermoType>
bool Foam::reactionModels::noReaction<ReactThermoType>::read()
{
    if (ReactThermoType::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
