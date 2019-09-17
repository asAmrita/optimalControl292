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
    scalar Jk = 0;
    scalar error = 0;

// Compute cost function value
#include "costFunctionValue.H"

    std::ofstream file("cost.csv");
    file << 0 << "," << J << "," << 0 << nl;

    std::ofstream errorFile("error.csv");

        runTime.write();

    while (runTime.loop() && fabs(J - Jold) > tol)
    {
        // save old cost value
        Jold = J;

        // Primal equation
        solve(fvm::laplacian(k, y)-fvm::Sp(1.0,y) + beta * u + f);
   
        // Adjoint equation
        solve(fvm::laplacian(k, p)-fvm::Sp(1.0,p) + y - yd);

        // Save current control
        uk = u;

        // calculate current cost
		#include "costFunctionValue.H"
        Jk = J;


        bool alphaFound = false;

        // calculate derivative^2 integrate((lambda*u + beta*p)^2 dv). Why??
        scalar phip0 = gSum(volField * Foam::pow(lambda * uk.internalField() + beta * p.internalField(), 2));

        while ((!alphaFound) && (alpha > tol))
        {
            u = uk - alpha * (lambda * uk + beta * p);

            // truncate u for constrained control set
            forAll(u, i)
            {
                u[i] = min(u[i], uMax[i]);
                u[i] = max(u[i], uMin[i]);
            }

            u.correctBoundaryConditions();

            // get new y
            solve(fvm::laplacian(k, y) -fvm::Sp(1.0, y) + beta * u + f);

            // get new cost
			#include "costFunctionValue.H"

            if (J <= Jk - c1 * alpha * phip0)
            {
                Info << "alpha found, alpha = " << alpha << ", J = " << J << ", phip0" << phip0 << endl;

                alphaFound = true;
            }
            else
            {
                Info << "alpha NOT found, alpha = " << alpha << endl;
                alpha = c2 * alpha;
            }
            Info<<J<<endl;
        }

        Info << "Iteration no. " << runTime.timeName() << " - "
             << "Cost value " << J
             << " - "
             << "Cost variation" << fabs(J - Jold) << endl;

        file.open("cost.csv",std::ios::app);
        file << runTime.value() << "," << J << nl;
        file.close();

        // calculate the L2 norm of the error in control variables
        error = Foam::sqrt(gSum(volField * (Foam::pow(u.internalField() - ud.internalField(), 2))));
        
        errorFile.open("error.csv",std::ios::app);
        errorFile << runTime.value() << "," << error << nl;
        errorFile.close();

        Info << "Iteration no. " << runTime.timeName() << " - " << "error " << error << endl;


       /* uDiff = u - ud;
        forAll(uDiff,i)
        {
            uDiff[i] = fabs(u[i] - ud[i]);
        }
      */
    uDiff=mag(u-ud);
    pDiff=mag(p-pd);
    yDiff=mag(y-yd);
        runTime.write();
    }

    file.close();
    errorFile.close();

    runTime++;
    y.write();
    p.write();
    u.write();
    uDiff.write();
    ud.write();
    //    uc.write();
    //    udiff.write();

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
