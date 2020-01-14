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

#include "chesubsurfaceTopFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chesubsurfaceTopFvPatchScalarField
::chesubsurfaceTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    YiFromSurf_(p.size(), 0.0),
    phiName_("phiSub")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::chesubsurfaceTopFvPatchScalarField
::chesubsurfaceTopFvPatchScalarField
(
    const chesubsurfaceTopFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    YiFromSurf_(ptf.YiFromSurf_, mapper),
    phiName_(ptf.phiName_)
{}


Foam::chesubsurfaceTopFvPatchScalarField
::chesubsurfaceTopFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    YiFromSurf_("YiFromSurf", dict, p.size()),
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


Foam::chesubsurfaceTopFvPatchScalarField
    ::chesubsurfaceTopFvPatchScalarField
(
    const chesubsurfaceTopFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    YiFromSurf_(ptf.YiFromSurf_),
    phiName_(ptf.phiName_)
{}


Foam::chesubsurfaceTopFvPatchScalarField
    ::chesubsurfaceTopFvPatchScalarField
(
    const chesubsurfaceTopFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    YiFromSurf_(ptf.YiFromSurf_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chesubsurfaceTopFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //fluid flux on this patch
    const fvsPatchField<scalar>& phiSubp =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    forAll(YiFromSurf_, facei)
    {
        if(phiSubp[facei] <= 0.0) //flux surface -> subsurface (inlet): fixedValue
        {
           this->refValue()[facei] = YiFromSurf_[facei];
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


void Foam::chesubsurfaceTopFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchField<scalar>::write(os);
    YiFromSurf_.writeEntry("YiFromSurf",os);
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
        chesubsurfaceTopFvPatchScalarField
    );
}

// ************************************************************************* //
