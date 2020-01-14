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
    surfaceBottomFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceBottomFixedValueFvPatchScalarField_H
#define surfaceBottomFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class surfaceBottomFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceBottomFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data
    // CFromSub: interpolated value of the same field from subsurface domain
    scalarField CFromSub_;
   
    // DsubFromSub_: interpolated value of dispersion coefficient from subsurface domain
    //scalarField DsubFromSub_;

    // deltaCoeffFromSub_: interploated value of delta (distance from cell center to boundary face center)
    // from subsurface domain
    //scalarField deltaCoeffFromSub_;

public:

    //- Runtime scalar information
    TypeName("surfaceBottomFixedValue");


    // Constructors

        //- Construct from patch and internal field
        surfaceBottomFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceBottomFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceBottomFixedValueFvPatchScalarField onto a new patch
        surfaceBottomFixedValueFvPatchScalarField
        (
            const surfaceBottomFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceBottomFixedValueFvPatchScalarField
        (
            const surfaceBottomFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new surfaceBottomFixedValueFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceBottomFixedValueFvPatchScalarField
        (
            const surfaceBottomFixedValueFvPatchScalarField&,
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
                new surfaceBottomFixedValueFvPatchScalarField(*this, iF)
            );
        }

    // Member functions
        virtual  scalarField& CFromSub()
        {
                return CFromSub_;
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
