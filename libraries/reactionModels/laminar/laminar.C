/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "laminar.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::reactionModels::laminar<Type>::laminar
(
    const word& modelType,
    const fvMesh& mesh,
    const word& reactionProperties,
    const word& phaseName
)
:
    Type(modelType, mesh, reactionProperties, phaseName),
    integrateReactionRate_
    (
        this->coeffs().lookupOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::reactionModels::laminar<Type>::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::reactionModels::laminar<Type>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class Type>
void Foam::reactionModels::laminar<Type>::correct()
{
    if (this->active())
    {
        if (integrateReactionRate_)
        {
            if (fv::localEulerDdt::enabled(this->mesh()))
            {
                const scalarField& rDeltaT =
                    fv::localEulerDdt::localRDeltaT(this->mesh());

                if (this->coeffs().found("maxIntegrationTime"))
                {
                    scalar maxIntegrationTime
                    (
                        readScalar(this->coeffs().lookup("maxIntegrationTime"))
                    );

                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxIntegrationTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvScalarMatrix>
Foam::reactionModels::laminar<Type>::R(volScalarField& Y) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimless*dimVolume/dimTime));

    fvScalarMatrix& Su = tSu.ref();

    if (this->active())
    {
        const label specieI =
            this->thermo().composition().species()[Y.member()];

        //Note: the += operator will multiply the RHS by cell volume
        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class Type>
bool Foam::reactionModels::laminar<Type>::read()
{
    if (Type::read())
    {
        integrateReactionRate_ =
            this->coeffs().lookupOrDefault("integrateReactionRate", true);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
