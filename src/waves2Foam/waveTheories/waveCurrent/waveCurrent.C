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

#include <functional> // For std::function
#include "waveCurrent.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveCurrent, 0);
addToRunTimeSelectionTable(waveTheory, waveCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveCurrent::waveCurrent(
    const word &subDictName,
    const fvMesh &mesh_)
    : waveTheory(subDictName, mesh_),
      H_(readScalar(coeffDict_.lookup("height"))),
      h_(readScalar(coeffDict_.lookup("depth"))),
      omega_(readScalar(coeffDict_.lookup("omega"))),
      period_(2 * PI_ / omega_),
      phi_(readScalar(coeffDict_.lookup("phi"))),
      k_(vector(coeffDict_.lookup("waveNumber"))),
      U_(vector(coeffDict_.lookup("U"))), // In uniform current, U is mean current velocity; in linear current, U is surface velocity; in log shear current, U is friction velocity.
    //   omegac_(omega_ + (k_ & U_)),
      K_(mag(k_)),
      Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", period_)),
      currentType_(coeffDict_.lookupOrDefault<word>("currentType", "uniform")), // valid current type: uniform, linear, log 
      nu_(coeffDict_.lookupOrDefault("nu", 1.002e-6)),
      ks_(coeffDict_.lookupOrDefault("ks", 1.25e-3)),
      kappa_(coeffDict_.lookupOrDefault<scalar>("kappa", 0.41)),
      Ad_(coeffDict_.lookupOrDefault<scalar>("Ad", 25.0)),
      wallDist_(wallDist::New(mesh_).y()), // distance with respect to bottom wall
      cellC_(mesh_.C()),
      debug_(Switch(coeffDict_.lookup("debug")))
    {
        checkWaveDirection(k_);

    if
    (
        H_/2.0 - 4.0*1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
        - 1.0/Foam::tanh(K_*h_)) < 0
    )
    {
        if (debug_)
        {
            WarningIn
            (
                "label waveCurrent::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
            Info << "a1 = " << H_/2.0 << " , a2 = "
                 << (1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                    - 1.0/Foam::tanh(K_*h_))) << endl;
        }
    }
    if (currentType_ == "uniform") {
        k1_ = waveNumberOne(h_, period_, U_.x(), U_.x());
    } else if (currentType_ == "linear") {
        k1_ = waveNumberOne(h_, period_, U_.x(), 0);
    } else if (currentType_ == "log") {
        k1_ = k_;
    } else {
        // Raise error, invalid current type.
        FatalErrorInFunction
            << "Invalid current type: " << currentType_ << ".\n"
            << "Valid current type: {uniform, linear, log}.\n" << exit(FatalError);
    }
    K1_ = mag(k1_);
}


void waveCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar waveCurrent::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2 * PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar waveCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar arg(omega_*time - (k1_ & x) + phi_);

    scalar eta = H_ / 2.0 * Foam::cos(arg) // First order contribution.
                 * factor(time)  // Hot-starting.
                 + seaLevel_;    // Adding sea level.

    return eta;
}


scalar waveCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
	scalar res = 0;

	// Get arguments and local coordinate system
    scalar Z(returnZ(x));
    scalar arg(omega_*time - (k1_ & x) + phi_);

    // First order contribution
    res = rhoWater_*mag(g_)*H_/2.0*Foam::cosh(K1_*(Z + h_))
        /Foam::cosh(K1_*h_)*Foam::cos(arg);
    // Apply the ramping-factor
    res *= factor(time);
    res += referencePressure();

    return res;
}



vector waveCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar cel(omega_/K1_); // FIXME: wtf?
    scalar arg(omega_*time - (k1_ & x) + phi_);

    // First order contribution
    scalar Uhorz = (H_ / 2) * (omega_ - (k1_ & U_)) *
                   Foam::cosh(K1_*(Z + h_))/Foam::cosh(K1_*h_) *
                   Foam::cos(arg);

    // Current velocity in x dirction
    Uhorz += currentU(x, time).x();

    // First order contribution
    // Note "-" because of "g" working in the opposite direction
    scalar Uvert = - (H_ / 2) * (omega_ - (k1_ & U_)) *
                   Foam::sinh(K1_*(Z + h_))/Foam::cosh(K1_*h_) *
                   Foam::sin(arg);

    // Current velocity in y dirction
    Uvert += currentU(x, time).y();

    // Multiply by the time stepping factor
    Uvert *= factor(time);
    Uhorz *= factor(time);

    // Generate the velocity vector
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k1_/K1_ - Uvert*direction_;
}

vector waveCurrent::currentU
(
    const point& x,
    const scalar& time
) const
{
    if (currentType_ == "uniform")
    {
        return U_;
    } else if (currentType_ == "linear")
    {
        // linear distributed current velocity, bottom is zero, top is 2U_
        volScalarField temp1=wallDist_;
        volScalarField temp2=wallDist_;

        temp1 *=scalar(0.0);
        temp2 *=scalar(0.0);
        scalar yPlus_ = 0.0;
        scalar yPlus_minus_ = 0.0;
        scalar distriFactor_ = 0.0; // distribution factor/weight
        vector shearU_ = vector(0,0,0);

        // Slow but consistent code
        forAll(temp1, index)
        {
            yPlus_ = wallDist_[index] / h_; // Not use, just for consistant

            distriFactor_ = 1.0;

            // https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2021/KorayDenizGoral/Report_KorayDenizGoral.pdf
            // page 47 line 147 code is wrong, 
            // temp1[0] or say temp[0] should be zero instead of vanDriest_(approx 1).
            temp2[index] = distriFactor_;

            if (index > 0) // Skip the first temp1 entry
            {
                yPlus_minus_ = wallDist_[index-1] / h_; // Not use, just for consistant
                temp1[index] = temp1[index-1] + 0.5*(temp2[index] + temp2[index-1])*(yPlus_ - yPlus_minus_);
            }

            if (cellC_[index] == x)
            {
                shearU_ = U_*temp1[index];
                return shearU_;
            }
        }
        // Never reach here
        FatalErrorInFunction
            << "Bad mesh, can not find cell at x = " << x << "." << exit(FatalError);
    } else if (currentType_ == "log")
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

        // Slow but consistent code
        forAll(temp1, index)
        {
            yPlus_ = wallDist_[index] * U_.x() / nu_;

            distriFactor_ = 1.0/(1.0+Foam::sqrt(1.0+4.0*Foam::pow(kappa_,2.0)*Foam::pow(yPlus_+deltayPlus_,2.0) \ 
                        *Foam::pow(1.0-Foam::exp(-(yPlus_+deltayPlus_)/Ad_),2.0)));

            // https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2021/KorayDenizGoral/Report_KorayDenizGoral.pdf
            // page 47 line 147 code is wrong, 
            // temp1[0] or say temp[0] should be zero instead of vanDriest_(approx 1).
            temp2[index] = distriFactor_;

            if (index > 0) // Skip the first temp1 entry
            {
                yPlus_minus_ = wallDist_[index-1] * U_.x() / nu_;
                temp1[index] = temp1[index-1] + 0.5*(temp2[index] + temp2[index-1])*(yPlus_ - yPlus_minus_);
            }

            if (cellC_[index] == x)
            {
                shearU_ = 2.0*U_*temp1[index];
                return shearU_;
            }
        }
        // Never reach here
        FatalErrorInFunction
            << "Bad mesh, can not find current velocity about cell at x = " << x << "." << exit(FatalError);
    } else 
    {
        // Raise error, invalid current type.
        FatalErrorInFunction
            << "Invalid current type: " << currentType_ << ".\n"
            << "Valid current type: {uniform, linear, log}.\n" << exit(FatalError);
    }
}

// Helper function for solving wavenumber one
scalar waveCurrent::secant_method(std::function<scalar(scalar)> func, scalar a, scalar b) const 
{
    double tol = 1e-6;
    int max_iter = 1000;
    scalar fa = func(a);
    scalar fb = func(b);
    for (int i = 0; i < max_iter; ++i) {
        if (abs(fa - fb) < 1e-10) { // Avoid division by zero
            return (a + b) / 2.0; // Approximate
        }
        scalar c = (a * fb - b * fa) / (fb - fa);
        scalar fc = func(c);
        if (abs(fc) < tol) {
            return c;
        }
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    // If no convergence, return approximation
    return (a + b) / 2.0;
}

vector waveCurrent::waveNumberOne
(
    const scalar h, // water depth
    const scalar T, // wave period
    const scalar Us, // surface velocity in x direction
    const scalar Ub // bottom velocity in x direction
) const
{
    scalar G = mag(g_);
    scalar k0 = k_.x(); // Initial guess of wave number
    scalar vorticity = (Us - Ub) / h;

    // Dispersion equation (lambda function)
    auto dispersion_eq = [=](scalar k) -> scalar {
        return pow(omega_ - k * Us, 2) -
               (G * k - vorticity * (omega_ - k * Us)) * tanh(k * h);
    };

    // Find solution in range [a, b]
    scalar a = k0 * 0.1;
    scalar b = k0 * 10.0;
    // Expand search range if not cross zero point
    if (dispersion_eq(a) * dispersion_eq(b) > 0) {
        if (dispersion_eq(k0) < 0) {
            a = k0 * 0.01;
        } else {
            b = k0 * 100.0;
        }
    }

    k1x = secant_method(dispersion_eq, a, b)
    return vector(k1x, 0, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
