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

#include "logCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace waveTheories
    {

        // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

        defineTypeNameAndDebug(logCurrent, 0);
        addToRunTimeSelectionTable(waveTheory, logCurrent, dictionary);

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


        logCurrent::logCurrent
            (
             const word& subDictName,
             const fvMesh& mesh_
            )
            :
                waveTheory(subDictName, mesh_),
                U_(vector(coeffDict_.lookup("U"))), // In uniform current, U is mean current velocity; in linear current, U is surface velocity; in log shear current, U is friction velocity.

                h_(readScalar(coeffDict_.lookup("depth"))),
                Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),
                nu_(coeffDict_.lookupOrDefault("nu", 1.002e-6)),
                ks_(coeffDict_.lookupOrDefault("ks", 1.25e-3)),
                kappa_(coeffDict_.lookupOrDefault<scalar>("kappa", 0.41)),
                Ad_(coeffDict_.lookupOrDefault<scalar>("Ad", 25.0)),
                cellC_(mesh_.C()), 
                wallDist_(wallDist::New(mesh_).y()), // distance with respect to bottom wall
                localSeaLevel_
                    (
                     coeffDict_.lookupOrDefault<scalar>("localSeaLevel", seaLevel_)
                    )
                    {}


        void logCurrent::printCoeffs()
        {
            Info << "Loading wave theory: " << typeName << endl;
        }


        // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


        scalar logCurrent::factor(const scalar& time) const
        {
            scalar factor(1);
            if (Tsoft_ > 0.0)
            {
                factor = Foam::sin(PI_/2.0/Tsoft_*Foam::min(Tsoft_, time));
            }

            return factor;
        }


        scalar logCurrent::eta
            (
             const point& x,
             const scalar& time
            ) const
            {
                //    scalar eta = seaLevel_;
                scalar eta = localSeaLevel_;
                return eta;
            }


        //scalar logCurrent::ddxPd
        //(
        //    const point& x,
        //    const scalar& time,
        //    const vector& unitVector
        //) const
        //{
        //    return 0.0;
        //}


        scalar logCurrent::pExcess
            (
             const point& x,
             const scalar& time
            ) const
            {
                return referencePressure(localSeaLevel_);
            }


        vector logCurrent::U
            (
             const point& x,
             const scalar& time
            ) const
            {
                // linear distributed current velocity, bottom is zero, top is 2U_
                volScalarField temp1=wallDist_;
                volScalarField temp2=wallDist_;

                temp1 *=scalar(0.0);
                temp2 *=scalar(0.0);
                scalar yPlus_ = 0.0;
                scalar yPlus_minus_ = 0.0;
                scalar distriFactor_ = 0.0; // distribution factor/weight

                scalar ksPlus_= ks_* U_.x() / nu_;
                scalar deltayPlus_ = 0.9*(Foam::sqrt(ksPlus_)-ksPlus_*Foam::exp(-ksPlus_/6.0));
                vector shearU_ = vector(0,0,0);

                forAll(temp1, index)
                {
                    yPlus_ = wallDist_[index] * U_.x() / nu_;

                    distriFactor_ = 1.0/(1.0+Foam::sqrt(1.0+4.0*Foam::pow(kappa_,2.0)*Foam::pow(yPlus_+deltayPlus_,2.0) \ 
                                *Foam::pow(1.0-Foam::exp(-(yPlus_+deltayPlus_)/Ad_),2.0))); // vanDriest_ factor

                    temp1[index] = distriFactor_; // near bed vanDriest_(approx 1).
                    temp2[index] = distriFactor_;

                    if (index > 0) // Skip the first temp1 entry
                    {
                        yPlus_minus_ = wallDist_[index-1] * U_.x() / nu_;
                        temp1[index] = temp1[index-1] + 0.5*(temp2[index] + temp2[index-1])*(yPlus_ - yPlus_minus_);
                    }

                    if (cellC_[index] == x) // Ignore the cells with 0 <= alpha < 1
                    {
                        shearU_ = 2.0*U_*temp1[index]; // Only alpha = 1 cells have non-zero current velocity
                        break;
                    }
                }
                return shearU_*factor(time);
            }


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
