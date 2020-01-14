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
    Foam::reactionModel

Description
    Base class for reaction models

SourceFiles
    reactionModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactionModel_H
#define reactionModel_H

#include "IOdictionary.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class reactionModel Declaration
\*---------------------------------------------------------------------------*/

class reactionModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Disallow copy construct
        reactionModel(const reactionModel&);

        //- Disallow default bitwise assignment
        void operator=(const reactionModel&);


protected:

    // Protected data

        //- Reference to the turbulence model
        incompressibleTurbulenceModel* turbulencePtr_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Active
        Switch active_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;

        //- Phase name
        const word phaseName_;


public:

    //- Runtime type information
    TypeName("reactionModel");

    //- Default reactionProperties dictionary name
    static const word reactionPropertiesName;


    // Constructors

        //- Construct from components
        reactionModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const word& reactionProperties=reactionPropertiesName,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~reactionModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Set turbulence
        inline void setTurbulence(incompressibleTurbulenceModel& turbModel);

        //- Return access to turbulence
        inline const incompressibleTurbulenceModel& turbulence() const;

        //- Return const access to rho
        inline const volScalarField& rho() const;

        //- Return const access to phi
        inline tmp<surfaceScalarField> phi() const;

        //- Is reaction active?
        inline const Switch& active() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct reaction rate
        virtual void correct() = 0;

        //- Specie consumption rate matrix, i.e. source term for specie equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
