/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceCO2FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceCO2FvPatchScalarField::surfaceCO2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    CO2s_(p.size(), 0.0),
    KCO2_L_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::surfaceCO2FvPatchScalarField::surfaceCO2FvPatchScalarField
(
    const surfaceCO2FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    CO2s_(ptf.CO2s_, mapper),
    KCO2_L_(ptf.KCO2_L_, mapper)
{}


Foam::surfaceCO2FvPatchScalarField::surfaceCO2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    CO2s_("CO2s", dict, p.size()),
    KCO2_L_("KCO2_L", dict, p.size())
{
    refValue() = CO2s_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::surfaceCO2FvPatchScalarField::surfaceCO2FvPatchScalarField
(
    const surfaceCO2FvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    CO2s_(tppsf.CO2s_),
    KCO2_L_(tppsf.KCO2_L_)
{}


Foam::surfaceCO2FvPatchScalarField::surfaceCO2FvPatchScalarField
(
    const surfaceCO2FvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    CO2s_(tppsf.CO2s_),
    KCO2_L_(tppsf.KCO2_L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceCO2FvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    CO2s_.autoMap(m);
    KCO2_L_.autoMap(m);
}


void Foam::surfaceCO2FvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const surfaceCO2FvPatchScalarField& tiptf =
        refCast<const surfaceCO2FvPatchScalarField>(ptf);

    CO2s_.rmap(tiptf.CO2s_, addr);
    KCO2_L_.rmap(tiptf.KCO2_L_, addr);
}


void Foam::surfaceCO2FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const label patchi = patch().index();

    valueFraction() =
        1.0/
        (
            1.0
          + turbModel.nuEff(patchi)*patch().deltaCoeffs()/KCO2_L_
        );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::surfaceCO2FvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    CO2s_.writeEntry("CO2s", os);
    KCO2_L_.writeEntry("KCO2_L", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        surfaceCO2FvPatchScalarField
    );
}

// ************************************************************************* //
