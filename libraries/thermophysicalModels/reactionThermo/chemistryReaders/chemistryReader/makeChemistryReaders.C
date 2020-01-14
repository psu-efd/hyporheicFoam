/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "makeReactionThermo.H"
#include "thermoPhysicsTypes.H"

#include "chemistryReader.H"
#include "foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// chemistry readers based on sensibleEnthalpy

//makeChemistryReader(constGasHThermoPhysics);
makeChemistryReader(constGenericHThermoPhysics);
//makeChemistryReader(constIncompressibleGasHThermoPhysics);

//makeChemistryReaderType(foamChemistryReader, constGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constGenericHThermoPhysics);
/*
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasHThermoPhysics
);
*/

// chemistry readers based on sensibleInternalEnergy

//makeChemistryReader(constGasEThermoPhysics);
makeChemistryReader(constGenericEThermoPhysics);
//makeChemistryReader(constIncompressibleGasEThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGenericEThermoPhysics);
/*
makeChemistryReaderType(foamChemistryReader, constGasEThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasEThermoPhysics
);
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

