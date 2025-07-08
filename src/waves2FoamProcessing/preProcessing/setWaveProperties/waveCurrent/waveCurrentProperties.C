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

#include "waveCurrentProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveCurrentProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    waveCurrentProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveCurrentProperties::waveCurrentProperties
(
    const Time& rT,
    dictionary& dict,
    vector g,
    bool write
)
:
    setWaveProperties(rT, dict, g, write),
    sfp_( rT, dict, g, false, "")
{
    Info << "\nConstructing: " << this->type() << endl;

    period_ = readScalar( dict.lookup("period") );
    depth_  = readScalar( dict.lookup("depth") );
    U_ = vector( dict_.lookup("U") ); // dict_ is same as dict, dict_ inherent from setWaveProperties, dict is passing argument.
    omega_  = 2.0*PI_/period_ ;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void waveCurrentProperties::set(Ostream& os)
{
    scalar k = sfp_.linearWaveNumber();
    scalar k1 = linearWaveNumberOne();

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );
    writeGiven( os, "U");

    if (dict_.found( "Tsoft" ))
    {
        writeGiven( os, "Tsoft");
    }

    if (dict_.found( "currentType" ))
    {
        writeGiven( os, "currentType");
    }

    writeGiven( os, "depth" );
    writeGiven( os, "period" );
    writeGiven( os, "direction" );
    writeGiven( os, "phi");
    writeGiven( os, "height");

    if (write_)
    {
        vector direction( vector(dict_.lookup("direction")));
        direction /= Foam::mag(direction);
        // direction *= k;

        writeDerived(os, "waveNumber", direction * k);
        writeDerived(os, "waveNumberOne", direction * k1);
        writeDerived(os, "omega", sfp_.omega());
    }

    writeGiven( os, "debug");

    // Write the relaxation zone
    writeRelaxationZone( os );

    // Write the closing bracket
    writeEnding( os );

    scalar H = readScalar( dict_.lookup("height") );
    scalar h = readScalar( dict_.lookup("depth")  );

    scalar a1 = H/2.0;
    scalar a2 = 1.0/16.0*k * sqr(H)
        *(3.0/Foam::pow(Foam::tanh(k*h),3.0) - 1.0/Foam::tanh(k*h));

    if (Switch( dict_.lookup("debug") ))
    {
        Info << nl << "The wave amplitudes are:\n" << tab << "  a1 = "
             << tab << a1
             << nl << tab << "  a2 = " << tab << a2
             << nl << tab << "4 a2 = " << tab << 4.0*a2
             << " (Validity criterion) " << endl;
    }

    if (a1 < 4.0*a2)
    {
        Info << a1 << tab << 4.0*a2 << endl;

        WarningIn
        (
            "void waveCurrentProperties::set(Ostream& os)"
        ) << endl << "The validity of Stokes second order is violated." << endl
          << endl;
    }
}

scalar waveCurrentProperties::linearWaveNumberOne() const
{
    scalar lower(0.0);

    scalar upper = Foam::max
        (
            4.0*PI_/( period_*Foam::sqrt( Foam::mag(G_)*depth_)),
            2.0*PI_/( Foam::pow( period_, 2.0))
        );

    scalar middle(0.5*(lower + upper) );

    scalar tanhMax(100);

    scalar valLower = Foam::pow(omega_ - middle * U_.x(), 2.0)
        - Foam::mag(G_)*lower*Foam::tanh( Foam::min(lower*depth_, tanhMax));
    scalar valUpper = Foam::pow(omega_ - middle * U_.x(), 2.0)
        - Foam::mag(G_)*upper*Foam::tanh( Foam::min(upper*depth_, tanhMax));
    scalar valMiddle = Foam::pow(omega_ - middle * U_.x(), 2.0)
        - Foam::mag(G_)*middle*Foam::tanh( Foam::min(middle*depth_, tanhMax));

    while (true)
    {
        if (Foam::sign( valLower ) == Foam::sign( valMiddle ))
        {
            lower    = middle;
            valLower = valMiddle;
        }
        else
        {
            upper    = middle;
            valUpper = valMiddle;
        }

        middle = 0.5*( lower + upper );

        valMiddle = Foam::pow(omega_ - middle * U_.x(), 2.0)
            - Foam::mag(G_)*middle
            *Foam::tanh( Foam::min(middle*depth_, tanhMax) );

        if
        (
            Foam::mag(valMiddle) < 1.0e-13 ||
            Foam::mag(valLower - valUpper)/middle < 1.0e-13
        )
        {
            break;
        }
    }

    return middle;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
