/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Application
    laplaceAdjointFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

    // Disable solvers performance output
    lduMatrix::debug = 0;
    solverPerformance::debug = 0;

    // Cost function value
    scalar J = 0;
    scalar Jold = 0;

// Compute cost function value
#include "costFunctionValue.H"

    std::ofstream file("results.csv");
    file << 0 << "," << J << "," << 0 << nl;

    while (runTime.loop() && fabs(J - Jold) > tol)
    {
        // Primal equation
        solve(fvm::laplacian(k, y) + beta * u);

        // Adjoint equation
        solve(fvm::laplacian(k, p) + y - yd);
        uc = -(beta/lambda)*p;
        udiff = u - uc;


        // Update control
        u = u - alpha * (lambda * u + beta * p);

        Jold = J;
#include "costFunctionValue.H"

        Info << "Iteration no. " << runTime.timeName() << " - "
             << "Cost value " << J
             << " - "
             << "Cost variation" << fabs(J - Jold) << endl;

        file << runTime.value() << "," << J << nl;

        runTime.write();
    }

    file.close();

    runTime++;
    y.write();
    p.write();
    u.write();
    uc.write();
    udiff.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << endl;
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "\nEnd\n"
         << endl;
    return 0;
}

// ************************************************************************* //
