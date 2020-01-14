/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::zeroNetFluxVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity boundary condition which ensures the net flux is zero

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                zeroNetFluxVelocity;
        value               uniform (0 0 0);
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    zeroNetFluxVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroNetFluxVelocityFvPatchVectorField_H
#define zeroNetFluxVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class zeroNetFluxVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class zeroNetFluxVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

    // Private member functions

        //- Update the patch values
        void updateValues();


public:

   //- Runtime type information
   TypeName("zeroNetFluxVelocity");


   // Constructors

        //- Construct from patch and internal field
        zeroNetFluxVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroNetFluxVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  zeroNetFluxVelocityFvPatchVectorField
        //  onto a new patch
        zeroNetFluxVelocityFvPatchVectorField
        (
            const zeroNetFluxVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroNetFluxVelocityFvPatchVectorField
        (
            const zeroNetFluxVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new zeroNetFluxVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroNetFluxVelocityFvPatchVectorField
        (
            const zeroNetFluxVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new zeroNetFluxVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
