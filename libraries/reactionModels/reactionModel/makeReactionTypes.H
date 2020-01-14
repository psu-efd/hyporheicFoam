/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#ifndef makeReactionTypes_H
#define makeReactionTypes_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactionTypesThermo(ReactModel, ReactType, Thermo, Table)          \
                                                                               \
    typedef Foam::reactionModels::ReactModel                                  \
        <Foam::reactionModels::ReactType, Foam::Thermo>                       \
        ReactModel##ReactType##Thermo;                                           \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        ReactModel##ReactType##Thermo,                                           \
        #ReactModel"<"#ReactType","#Thermo">",                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace reactionModels                                             \
        {                                                                      \
            typedef ReactModel<ReactType, Thermo> ReactModel##ReactType##Thermo;   \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                ReactModel##ReactType##Thermo,                                   \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


#define makeReactionTypes(ReactModel, ReactType, Table)                        \
                                                                               \
    typedef Foam::reactionModels::ReactModel                                  \
        <Foam::reactionModels::ReactType>                                     \
        ReactModel##ReactType;                                                   \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        ReactModel##ReactType,                                                   \
        #ReactModel"<"#ReactType">",                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace reactionModels                                             \
        {                                                                      \
            typedef ReactModel<ReactType> ReactModel##ReactType;                   \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Table,                                                         \
                ReactModel##ReactType,                                           \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
