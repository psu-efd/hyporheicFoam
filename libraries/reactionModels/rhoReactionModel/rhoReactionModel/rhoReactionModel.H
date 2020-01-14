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
    Foam::rhoReactionModel

Description
    Reaction models for rho-based thermodynamics

SourceFiles
    rhoReactionModelI.H
    rhoReactionModel.C
    rhoReactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef rhoReactionModel_H
#define rhoReactionModel_H

#include "reactionModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                     class rhoReactionModel Declaration
\*---------------------------------------------------------------------------*/

class rhoReactionModel
:
    public reactionModel
{

    // Private Member Functions

        //- Construct as copy (not implemented)
        rhoReactionModel(const rhoReactionModel&);

        //- Disallow default bitwise assignment
        void operator=(const rhoReactionModel&);


public:

    typedef rhoReactionThermo ReactionThermo;


    //- Runtime type information
    TypeName("rhoReactionModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoReactionModel,
        dictionary,
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& reactionProperties,
            const word& phaseName
        ),
        (modelType, mesh, reactionProperties, phaseName)
    );


    // Constructors


        //- Construct from components
        rhoReactionModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& reactionProperties,
            const word& phaseName
        );



    //- Selector
    static autoPtr<rhoReactionModel> New
    (
        const fvMesh& mesh,
        const word& reactionProperties=reactionPropertiesName,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~rhoReactionModel();


    // Member Functions

        // Access functions

            //- Access reaction dict
             inline const dictionary& coeff() const;


        //- Return access to the thermo package
        virtual rhoReactionThermo& thermo() = 0;

        //- Return const access to the thermo package
        virtual const rhoReactionThermo& thermo() const = 0;


     // IO

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
