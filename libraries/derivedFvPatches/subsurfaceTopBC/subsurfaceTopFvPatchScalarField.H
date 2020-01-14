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
    Foam::subsurfaceTopFvPatchScalarField


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
    subsurfaceTopFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef subsurfaceTopFvPatchScalarField_H
#define subsurfaceTopFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class subsurfaceTopFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class subsurfaceTopFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data
    // CFromSurf: interpolated value of the same field from surface domain
    scalarField CFromSurf_;

    //- Name of flux field
    word phiName_;

 
public:

    //- Runtime scalar information
    TypeName("subsurfaceTop");


    // Constructors

        //- Construct from patch and internal field
        subsurfaceTopFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subsurfaceTopFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  subsurfaceTopFvPatchScalarField onto a new patch
        subsurfaceTopFvPatchScalarField
        (
            const subsurfaceTopFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subsurfaceTopFvPatchScalarField
        (
            const subsurfaceTopFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new subsurfaceTopFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subsurfaceTopFvPatchScalarField
        (
            const subsurfaceTopFvPatchScalarField&,
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
                new subsurfaceTopFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        virtual  scalarField& CFromSurf()
        {
                return CFromSurf_;
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
