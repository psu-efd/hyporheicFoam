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
    Foam::surfaceBottomFvPatchScalarField


Description
    A scalar BC condition for the bottom patch of the surface flow domain.
    
    This BC relies on the flux direction on the bottom patch:
    
          (1) if flux direction subsurface -> surface (inlet)
                      fixedValue: the value comes from the subsurface (CFromSub_)

          (2) if flux direction subsurface <- surface (outlet)
                      zeroGradient: just let the scalar go out

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            surfaceBottom;
        phi             phi;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    surfaceBottomFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceBottomFvPatchScalarField_H
#define surfaceBottomFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class surfaceBottomFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class surfaceBottomFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data
    // CFromSub: interpolated value of the same field from subsurface domain
    scalarField CFromSub_;

    //- Name of flux field
    word phiName_;


public:

    //- Runtime scalar information
    TypeName("surfaceBottom");


    // Constructors

        //- Construct from patch and internal field
        surfaceBottomFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceBottomFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceBottomFvPatchScalarField onto a new patch
        surfaceBottomFvPatchScalarField
        (
            const surfaceBottomFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        surfaceBottomFvPatchScalarField
        (
            const surfaceBottomFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new surfaceBottomFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        surfaceBottomFvPatchScalarField
        (
            const surfaceBottomFvPatchScalarField&,
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
                new surfaceBottomFvPatchScalarField(*this, iF)
            );
        }

    // Member functions
        virtual  scalarField& CFromSub()
        {
                return CFromSub_;
        }

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
