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

#include "zeroNetFluxVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroNetFluxVelocityFvPatchVectorField::
zeroNetFluxVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{}


Foam::zeroNetFluxVelocityFvPatchVectorField::
zeroNetFluxVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false)
{
    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::zeroNetFluxVelocityFvPatchVectorField::
zeroNetFluxVelocityFvPatchVectorField
(
    const zeroNetFluxVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::zeroNetFluxVelocityFvPatchVectorField::
zeroNetFluxVelocityFvPatchVectorField
(
    const zeroNetFluxVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{}


Foam::zeroNetFluxVelocityFvPatchVectorField::
zeroNetFluxVelocityFvPatchVectorField
(
    const zeroNetFluxVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroNetFluxVelocityFvPatchVectorField::updateValues()
{
    const vectorField n(patch().nf());

    vectorField Up(*this);

    // Patch normal velocity
    scalarField nUp(n & Up);

    // Remove the normal component of the extrapolate patch velocity
    Up -= nUp*n;

    const scalar flowRate = 0.0;
    const scalar estimatedFlowRate = -gSum((this->patch().magSf()*nUp));

    nUp -= ((flowRate - estimatedFlowRate)/gSum(patch().magSf()));

    // Add the corrected normal component of velocity to the patch velocity
    Up += nUp*n;

    // Correct the patch velocity
    this->operator==(Up);

    //Info << "zeroNetFluxVelocity = " << Up << endl;
}


void Foam::zeroNetFluxVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateValues();

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::zeroNetFluxVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       zeroNetFluxVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
