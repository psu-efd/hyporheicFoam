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
    Foam::surfaceCO2FvPatchScalarField

Group
    grpWallBoundaryConditions

Description
    This boundary condition provides a boundary conditon for reaeration
    through the free surface

Usage
    \table
        Property     | Description              | Required    | Default value
        CO2s           | saturated CO2 conc.       | yes         |
        KCO2_L          | liquid-film coeff. [L/T] | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceCO2;
        CO2s              uniform 0.006;    // 1 kg/m^3 = 1000 g/m^3 = 1000 mg/L
        KCO2_L             uniform 3.47e-6;  // m/s  (=Ka*H)        
    }
    \endverbatim

SourceFiles
    surfaceCO2FvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceCO2FvPatchScalarField_H
#define surfaceCO2FvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class surfaceCO2FvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceCO2FvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Saturated CO2 concentration
        scalarField CO2s_;

        //- Liquid-film coeff.
        scalarField KCO2_L_;


public:

    //- Runtime type information
    TypeName("surfaceCO2");


    // Constructors

        //- Construct from patch and internal field
        surfaceCO2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceCO2FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given surfaceCO2FvPatchScalarField
        // onto a new patch
        surfaceCO2FvPatchScalarField
        (
            const surfaceCO2FvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceCO2FvPatchScalarField
        (
            const surfaceCO2FvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new surfaceCO2FvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceCO2FvPatchScalarField
        (
            const surfaceCO2FvPatchScalarField&,
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
                new surfaceCO2FvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return CO2s
            const scalarField& CO2s() const
            {
                return CO2s_;
            }

            //- Return reference to CO2s to allow adjustment
            scalarField& CO2s()
            {
                return CO2s_;
            }

            //- Return KCO2_L
            const scalarField& KCO2_L() const
            {
                return KCO2_L_;
            }

            //- Return reference to KCO2_L to allow adjustment
            scalarField& KCO2_L()
            {
                return KCO2_L_;
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
