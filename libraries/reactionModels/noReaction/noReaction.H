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
    Foam::reactionModels::noReaction

Description
    Dummy reaction model for 'no reaction'

SourceFiles
    noReaction.C

\*---------------------------------------------------------------------------*/

#ifndef noReaction_H
#define noReaction_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                        Class noReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactThermoType>
class noReaction
:
    public ReactThermoType
{

    //- Disallow copy construct
    noReaction(const noReaction&);

    //- Disallow default bitwise assignment
    void operator=(const noReaction&);


public:

    //- Runtime type information
    TypeName("noReaction");


    // Constructors

        //- Construct from components
        noReaction
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& reactionProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~noReaction();


    // Member Functions

        //- Correct reaction rate
        virtual void correct();

        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
