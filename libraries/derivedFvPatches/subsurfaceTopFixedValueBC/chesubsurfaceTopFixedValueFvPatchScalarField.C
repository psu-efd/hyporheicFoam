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

#include "chesubsurfaceTopFixedValueFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chesubsurfaceTopFixedValueFvPatchScalarField
::chesubsurfaceTopFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    YiFromSurf_(p.size(), 0.0)
    //DsubFromSub_(p.size(), 0.0)
   // deltaCoeffFromSub_(p.size(), 0.0)
   // realFlux_(p.size(), 0.0),
   // Ic_(p.size(), 0.0)
{
}


Foam::chesubsurfaceTopFixedValueFvPatchScalarField
::chesubsurfaceTopFixedValueFvPatchScalarField
(
    const chesubsurfaceTopFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    YiFromSurf_(ptf.YiFromSurf_, mapper)
    //DsubFromSub_(ptf.DsubFromSub_, mapper)
    //deltaCoeffFromSub_(ptf.deltaCoeffFromSub_, mapper)
 //   realFlux_(ptf.realFlux_, mapper),
 //   Ic_(ptf.Ic_, mapper)
{}


Foam::chesubsurfaceTopFixedValueFvPatchScalarField
::chesubsurfaceTopFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    YiFromSurf_("YiFromSurf", dict, p.size())
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


Foam::chesubsurfaceTopFixedValueFvPatchScalarField
    ::chesubsurfaceTopFixedValueFvPatchScalarField
(
    const chesubsurfaceTopFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    YiFromSurf_(ptf.YiFromSurf_)
    //DsubFromSub_(ptf.DsubFromSub_)
    //deltaCoeffFromSub_(ptf.deltaCoeffFromSub_)
//    realFlux_(ptf.realFlux_),
//    Ic_(ptf.Ic_)
{}


Foam::chesubsurfaceTopFixedValueFvPatchScalarField
    ::chesubsurfaceTopFixedValueFvPatchScalarField
(
    const chesubsurfaceTopFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    YiFromSurf_(ptf.YiFromSurf_)
    //DsubFromSub_(ptf.DsubFromSub_)
    //deltaCoeffFromSub_(ptf.deltaCoffFromSub_)
//    realFlux_(ptf.realFlux_),
//    Ic_(ptf.Ic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::chesubsurfaceTopFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //scalar value on the top patch of the subsurface doamin
//    const fvPatchField<scalar>& CFromSubp =
//        patch().lookupPatchField<volScalarField, scalar>("CFromSub");

    //dispersion coefficient on the top patch of the subsurface domain
    vectorField Sn_f = patch().Sf()/patch().magSf();

    //const fvPatchField<scalar>& DFromSubp =
    const scalarField DFromSubp =
        ((Sn_f & (patch().lookupPatchField<volTensorField, scalar>("DSub_Coeff"))) & Sn_f);

    //diffusion coefficient in the surface domain due to turbulence
    const fvPatchField<scalar>& DtFromSurfp =
        patch().lookupPatchField<volScalarField, scalar>("DFromSurf");

    //delta coefficient of the top patch in the subsurface domain
    const fvPatchField<scalar>& deltaCoeffFromSurfp =
        patch().lookupPatchField<volScalarField, scalar>("deltaCoeffFromSurf");

    const scalarField deltaCoeffFromSubp(this->patch().deltaCoeffs());

    //porosity
    scalar theta = 0.33;

/*
    Info << "DtFromSurfp = " << DtFromSurfp << endl;
    Info << "DFromSubp = " << DFromSubp << endl;
    Info << "deltaCoeffFromSurfp = " << deltaCoeffFromSurfp << endl;
    Info << "deltaCoeffFromSubp = " << deltaCoeffFromSubp << endl;
*/
    //note: here it should be multiplied by the porosity theta to get the 
    //bulk concentration.
    operator==
    (
       theta*(DtFromSurfp/deltaCoeffFromSubp*YiFromSurf_ + DFromSubp/deltaCoeffFromSurfp*(this->patchInternalField()))
      /(DtFromSurfp/deltaCoeffFromSubp + theta*DFromSubp/deltaCoeffFromSurfp)
    );

    //Info << "subsurface top value = " << (*this) << endl;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::chesubsurfaceTopFixedValueFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
    YiFromSurf_.writeEntry("YiFromSurf",os);
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
        chesubsurfaceTopFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //
