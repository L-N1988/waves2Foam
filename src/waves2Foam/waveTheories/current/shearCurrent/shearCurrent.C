/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "shearCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(shearCurrent, 0);
addToRunTimeSelectionTable(waveTheory, shearCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


shearCurrent::shearCurrent
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    U_(vector(coeffDict_.lookup("U"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),
    cellC_(mesh_.C()), 
    wallDist_(wallDist::New(mesh_).y()), // distance with respect to bottom wall
    localSeaLevel_
    (
        coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
    )
{}


void shearCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar shearCurrent::factor(const scalar& time) const
{
    scalar factor(1);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(PI_/2.0/Tsoft_*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar shearCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
//    scalar eta = seaLevel_;
    scalar eta = localSeaLevel_;
    return eta;
}


//scalar shearCurrent::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    return 0.0;
//}


scalar shearCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
    return referencePressure(localSeaLevel_);
}


vector shearCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    // linear distributed current velocity, bottom is zero, top is U_
    vector shearU_ = vector(0,0,0);

    forAll(cellC_, index)
    {
        if (cellC_[index] == x) // Only apply for submerged cells (alpha = 1), not for partially filled (0 < alpha < 1) or empty cells (alpha = 0)
        {
            // The partially filled cells can nerver reached here,
            // since the point x are re-calculated based on negative (wet) portion cells,
            // which can not find matched centroid of any cell in mesh.
            // The input point x comes from line 125 relaxationSchemeSpatial.C:
            // UTarget = waveProps_->U(lc.centreNeg(), mesh_.time().value());
            // lc means local cell, Neg means wet portion
            shearU_ = U_ * wallDist_[index] / h_;
            break;
        }
    }
    return shearU_*factor(time);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
