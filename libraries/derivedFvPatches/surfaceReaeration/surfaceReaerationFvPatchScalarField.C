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

#include "surfaceReaerationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceReaerationFvPatchScalarField::surfaceReaerationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Cs_(p.size(), 0.0),
    K_L_(p.size(), 0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::surfaceReaerationFvPatchScalarField::surfaceReaerationFvPatchScalarField
(
    const surfaceReaerationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Cs_(ptf.Cs_, mapper),
    K_L_(ptf.K_L_, mapper)
{}


Foam::surfaceReaerationFvPatchScalarField::surfaceReaerationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Cs_("Cs", dict, p.size()),
    K_L_("K_L", dict, p.size())
{
    refValue() = Cs_;
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


Foam::surfaceReaerationFvPatchScalarField::surfaceReaerationFvPatchScalarField
(
    const surfaceReaerationFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    Cs_(tppsf.Cs_),
    K_L_(tppsf.K_L_)
{}


Foam::surfaceReaerationFvPatchScalarField::surfaceReaerationFvPatchScalarField
(
    const surfaceReaerationFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    Cs_(tppsf.Cs_),
    K_L_(tppsf.K_L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceReaerationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Cs_.autoMap(m);
    K_L_.autoMap(m);
}


void Foam::surfaceReaerationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const surfaceReaerationFvPatchScalarField& tiptf =
        refCast<const surfaceReaerationFvPatchScalarField>(ptf);

    Cs_.rmap(tiptf.Cs_, addr);
    K_L_.rmap(tiptf.K_L_, addr);
}


void Foam::surfaceReaerationFvPatchScalarField::updateCoeffs()
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
          + turbModel.nuEff(patchi)*patch().deltaCoeffs()/K_L_
        );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::surfaceReaerationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    Cs_.writeEntry("Cs", os);
    K_L_.writeEntry("K_L", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        surfaceReaerationFvPatchScalarField
    );
}

// ************************************************************************* //
