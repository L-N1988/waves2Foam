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

  Class
  Foam::waveTheories::logCurrent

  Description
  Implementation of a shear current velocity profile, e.g. linear
  velocity distribution over the vertical. In the special case of
  U = vector::zero it can be used as outlet sponge layer wave type.

  A description of the general wave framework is given in

  @article { jacobsenFuhrmanFredsoe2011,
  Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
  title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
  Journal = {{Int. J. for Numer. Meth. Fluids}},
  Year = {2012},
  Volume = {70},
  Number = {9},
  Pages = {1073-1088},
  DOI = {{10.1002/fld.2726}},
  }

  SourceFiles
  logCurrent.C

  Author
  Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.

  \*---------------------------------------------------------------------------*/

#ifndef logCurrent_H
#define logCurrent_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace waveTheories
    {

        /*---------------------------------------------------------------------------*\
          Class logCurrent Declaration
          \*---------------------------------------------------------------------------*/

        class logCurrent
            :
                public waveTheory
        {

            protected:

                // Protected data
                vector U_;
                scalar h_; // water depth
                scalar Tsoft_;
                scalar nu_;
                scalar ks_;
                scalar kappa_;
                scalar Ad_;
                const volScalarField& wallDist_;
                const volVectorField& cellC_;
                scalar localSeaLevel_;

                // Protected member functions

                scalar factor(const scalar&) const;
            public:

                //- Runtime type information
                TypeName("logCurrent");

                // Constructors

                //- from components
                logCurrent
                    (
                     const word&,
                     const fvMesh& mesh_
                    );


                // Destructor

                ~logCurrent()
                {}


                // Member Functions
                virtual Switch conflictTSoftInitialise() const
                {
                    if (Foam::mag(U_) > SMALL && Tsoft_ > 0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                };

                void printCoeffs();

                virtual scalar eta
                    (
                     const point&,
                     const scalar&
                    ) const;

                //        virtual scalar ddxPd
                //        (
                //            const point&,
                //            const scalar&,
                //            const vector&
                //        ) const;

                virtual scalar pExcess
                    (
                     const point&,
                     const scalar&
                    ) const;

                bool implementPressure() const
                {
                    return true;
                };

                virtual vector U
                    (
                     const point&,
                     const scalar&
                    ) const;

        };


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
