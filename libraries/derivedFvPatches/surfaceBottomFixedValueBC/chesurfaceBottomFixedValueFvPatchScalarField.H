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

Class
    Foam::surfaceBottomFixedValueFvPatchScalarField


Description
    A scalar BC condition for the bottom patch of the surface flow domain.
    
    This BC uses a fixed value which is a weighted mean of the surface and subsurface cell values
    to ensure the same diffusive flux.
    
    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            surfaceBottom;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    chesurfaceBottomFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chesurfaceBottomFixedValueFvPatchScalarField_H
#define chesurfaceBottomFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class chesurfaceBottomFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chesurfaceBottomFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data
    // CFromSub: interpolated value of the same field from subsurface domain
    scalarField YiFromSub_;
   
    // DsubFromSub_: interpolated value of dispersion coefficient from subsurface domain
    //scalarField DsubFromSub_;

    // deltaCoeffFromSub_: interploated value of delta (distance from cell center to boundary face center)
    // from subsurface domain
    //scalarField deltaCoeffFromSub_;

public:

    //- Runtime scalar information
    TypeName("chesurfaceBottomFixedValue");


    // Constructors

        //- Construct from patch and internal field
        chesurfaceBottomFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chesurfaceBottomFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceBottomFixedValueFvPatchScalarField onto a new patch
        chesurfaceBottomFixedValueFvPatchScalarField
        (
            const chesurfaceBottomFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        chesurfaceBottomFixedValueFvPatchScalarField
        (
            const chesurfaceBottomFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new chesurfaceBottomFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        chesurfaceBottomFixedValueFvPatchScalarField
        (
            const chesurfaceBottomFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new chesurfaceBottomFixedValueFvPatchScalarField(*this, iF)
            );
        }

    // Member functions
        virtual  scalarField& YiFromSub()
        {
                return YiFromSub_;
        }
/*
        virtual  scalarField& DsubFromSub()
        {
                return DsubFromSub_;
        }

        virtual  scalarField& deltaCoeffFromSub()
        {
                return deltaCoeffFromSub_;
        }
*/

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
