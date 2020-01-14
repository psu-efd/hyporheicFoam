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

Class
    Foam::surfaceReaerationFvPatchScalarField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a boundary conditon for reaeration
    through the free surface

Usage
    \table
        Property     | Description              | Required    | Default value
        Cs           | saturated O2 conc.       | yes         |
        K_L          | liquid-film coeff. [L/T] | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceReaeration;
        Cs              uniform 0.006;    // 1 kg/m^3 = 1000 g/m^3 = 1000 mg/L
        K_L             uniform 3.47e-6;  // m/s  (=Ka*H)        
    }
    \endverbatim

SourceFiles
    surfaceReaerationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceReaerationFvPatchScalarField_H
#define surfaceReaerationFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class surfaceReaerationFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceReaerationFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Saturated O2 concentration
        scalarField Cs_;

        //- Liquid-film coeff.
        scalarField K_L_;


public:

    //- Runtime type information
    TypeName("surfaceReaeration");


    // Constructors

        //- Construct from patch and internal field
        surfaceReaerationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceReaerationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given surfaceReaerationFvPatchScalarField
        // onto a new patch
        surfaceReaerationFvPatchScalarField
        (
            const surfaceReaerationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceReaerationFvPatchScalarField
        (
            const surfaceReaerationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceReaerationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceReaerationFvPatchScalarField
        (
            const surfaceReaerationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceReaerationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return Cs
            const scalarField& Cs() const
            {
                return Cs_;
            }

            //- Return reference to Cs to allow adjustment
            scalarField& Cs()
            {
                return Cs_;
            }

            //- Return K_L
            const scalarField& K_L() const
            {
                return K_L_;
            }

            //- Return reference to K_L to allow adjustment
            scalarField& K_L()
            {
                return K_L_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

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
