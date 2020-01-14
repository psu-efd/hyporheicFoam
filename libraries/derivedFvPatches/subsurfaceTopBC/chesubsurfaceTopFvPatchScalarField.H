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
    Foam::chesubsurfaceTopFvPatchScalarField


Description
    A scalar BC condition for the top patch of the subsurface flow domain.
    
    This BC relies on the flux direction on subsurface top:
    
          (1) if flux direction subsurface -> surface (outlet)
                      zeroGradient: just let the scalar go out

          (2) if flux direction subsurface <- surface (inlet)
                      fixedValue: the value comes from the surface

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            subsurfaceTop;
        phi             phiSub;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    chesubsurfaceTopFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef chesubsurfaceTopFvPatchScalarField_H
#define chesubsurfaceTopFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class chesubsurfaceTopFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class chesubsurfaceTopFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data
    // CFromSurf: interpolated value of the same field from surface domain
    scalarField YiFromSurf_;

    //- Name of flux field
    word phiName_;

 
public:

    //- Runtime scalar information
    TypeName("chesubsurfaceTop");


    // Constructors

        //- Construct from patch and internal field
        chesubsurfaceTopFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        chesubsurfaceTopFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  subsurfaceTopFvPatchScalarField onto a new patch
        chesubsurfaceTopFvPatchScalarField
        (
            const chesubsurfaceTopFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        chesubsurfaceTopFvPatchScalarField
        (
            const chesubsurfaceTopFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new chesubsurfaceTopFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        chesubsurfaceTopFvPatchScalarField
        (
            const chesubsurfaceTopFvPatchScalarField&,
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
                new chesubsurfaceTopFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        virtual  scalarField& YiFromSurf()
        {
                return YiFromSurf_;
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
