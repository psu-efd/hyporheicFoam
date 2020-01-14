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

#include "chesurfaceBottomFixedValueFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chesurfaceBottomFixedValueFvPatchScalarField
::chesurfaceBottomFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    YiFromSub_(p.size(), 0.0)
    //DsubFromSub_(p.size(), 0.0)
   // deltaCoeffFromSub_(p.size(), 0.0)
   // realFlux_(p.size(), 0.0),
   // Ic_(p.size(), 0.0)
{
}


Foam::chesurfaceBottomFixedValueFvPatchScalarField
::chesurfaceBottomFixedValueFvPatchScalarField
(
    const chesurfaceBottomFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    YiFromSub_(ptf.YiFromSub_, mapper)
    //DsubFromSub_(ptf.DsubFromSub_, mapper)
    //deltaCoeffFromSub_(ptf.deltaCoeffFromSub_, mapper)
 //   realFlux_(ptf.realFlux_, mapper),
 //   Ic_(ptf.Ic_, mapper)
{}


Foam::chesurfaceBottomFixedValueFvPatchScalarField
::chesurfaceBottomFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    YiFromSub_("YiFromSub", dict, p.size())
    //DsubFromSub_("DsubFromSub", dict, p.size())
    //deltaCoeffFromSub_("deltaCoeffFromSub", dict, p.size())
 //   realFlux_("realFlux", dict, p.size()),
 //   Ic_("Ic", dict, p.size())
{

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
}


Foam::chesurfaceBottomFixedValueFvPatchScalarField
    ::chesurfaceBottomFixedValueFvPatchScalarField
(
    const chesurfaceBottomFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    YiFromSub_(ptf.YiFromSub_)
    //DsubFromSub_(ptf.DsubFromSub_)
    //deltaCoeffFromSub_(ptf.deltaCoeffFromSub_)
//    realFlux_(ptf.realFlux_),
//    Ic_(ptf.Ic_)
{}


Foam::chesurfaceBottomFixedValueFvPatchScalarField
    ::chesurfaceBottomFixedValueFvPatchScalarField
(
    const chesurfaceBottomFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    YiFromSub_(ptf.YiFromSub_)
    //DsubFromSub_(ptf.DsubFromSub_)
    //deltaCoeffFromSub_(ptf.deltaCoffFromSub_)
//    realFlux_(ptf.realFlux_),
//    Ic_(ptf.Ic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chesurfaceBottomFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //scalar value on the top patch of the subsurface doamin
//    const fvPatchField<scalar>& CFromSubp =
//        patch().lookupPatchField<volScalarField, scalar>("CFromSub");

    //dispersion coefficient on the top patch of the subsurface domain
    const fvPatchField<scalar>& DFromSubp =
        patch().lookupPatchField<volScalarField, scalar>("DFromSub");

    //diffusion coefficient in the surface domain due to turbulence
    const fvPatchField<scalar>& DtFromSurfp =
        patch().lookupPatchField<volScalarField, scalar>("Dt_Coeff");

    //delta coefficient of the top patch in the subsurface domain
    const fvPatchField<scalar>& deltaCoeffFromSubp =
        patch().lookupPatchField<volScalarField, scalar>("deltaCoeffFromSub");

    const scalarField deltaCoeffFromSurfp(this->patch().deltaCoeffs());

    //porosity
    scalar theta = 0.33;

/*
    Info << "deltaCoeffFromSubp = " << deltaCoeffFromSubp << endl;
    Info << "deltaCoeffFromSurfp = " << deltaCoeffFromSurfp << endl;
    Info << "DtFromSurfp = " << DtFromSurfp << endl;
    Info << "DFromSubp = " << DFromSubp << endl;
    Info << "CFromSurf = " << this->patchInternalField() << endl;
    Info << "CFromSub_ = " << CFromSub_ << endl;
*/
    operator==
    (
       (DtFromSurfp/deltaCoeffFromSubp*(this->patchInternalField()) + DFromSubp/deltaCoeffFromSurfp*YiFromSub_)
      /(DtFromSurfp/deltaCoeffFromSubp + theta*DFromSubp/deltaCoeffFromSurfp)
    );

    //Info << "surface bottom value = " << (*this) << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::chesurfaceBottomFixedValueFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    YiFromSub_.writeEntry("YiFromSub",os);
    //DsubFromSub_.writeEntry("DsubFromSub",os);
    //deltaCoeffFromSub_.writeEntry("deltaCoeffFromSub",os);
   // realFlux_.writeEntry("realFlux",os);
   // Ic_.writeEntry("Ic",os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        chesurfaceBottomFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //
