/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "subsurfaceTopFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subsurfaceTopFvPatchScalarField
::subsurfaceTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    CFromSurf_(p.size(), 0.0),
    phiName_("phiSub")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::subsurfaceTopFvPatchScalarField
::subsurfaceTopFvPatchScalarField
(
    const subsurfaceTopFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    CFromSurf_(ptf.CFromSurf_, mapper),
    phiName_(ptf.phiName_)
{}


Foam::subsurfaceTopFvPatchScalarField
::subsurfaceTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    CFromSurf_("CFromSurf", dict, p.size()),
    phiName_(dict.lookupOrDefault<word>("phi", "phiSub"))
{
    this->refValue() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->patchInternalField());
    }

    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::subsurfaceTopFvPatchScalarField
    ::subsurfaceTopFvPatchScalarField
(
    const subsurfaceTopFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    CFromSurf_(ptf.CFromSurf_),
    phiName_(ptf.phiName_)
{}


Foam::subsurfaceTopFvPatchScalarField
    ::subsurfaceTopFvPatchScalarField
(
    const subsurfaceTopFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    CFromSurf_(ptf.CFromSurf_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subsurfaceTopFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //fluid flux on this patch
    const fvsPatchField<scalar>& phiSubp =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    forAll(CFromSurf_, facei)
    {
        if(phiSubp[facei] <= 0.0) //flux surface -> subsurface (inlet): fixedValue
        {
           this->refValue()[facei] = CFromSurf_[facei];
           this->valueFraction()[facei] = 1.0;
        }
        else                      //flux surface <- subsurface (outlet): zeroGradient
        {
           this->refGrad()[facei] = 0.0;
           this->valueFraction()[facei] = 0.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::subsurfaceTopFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    CFromSurf_.writeEntry("CFromSurf",os);
    if (phiName_ != "phiSub")
    {
        os.writeKeyword("phiSub") << phiName_ << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        subsurfaceTopFvPatchScalarField
    );
}

// ************************************************************************* //
