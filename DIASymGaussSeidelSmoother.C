/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(a)
#define LIKWID_MARKER_STOP(a)
#define LIKWID_MARKER_CLOSE
#endif

#ifdef LIKWID_PERFMON
__attribute__((constructor)) static void likwid_init_wrapper()
{
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_REGISTER("DIASYMsweep");
}

__attribute__((destructor)) static void likwid_close_wrapper()
{
    LIKWID_MARKER_CLOSE;
}
#endif

#include "DIASymGaussSeidelSmoother.H"
#include <set>
#include <vector>
#include "symGaussSeidelSmoother.H"
#include <chrono>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DIASymGaussSeidelSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<DIASymGaussSeidelSmoother>
        addDIASymGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<DIASymGaussSeidelSmoother>
        addDIASymGaussSeidelSmootherAsymMatrixConstructorToTable_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DIASymGaussSeidelSmoother::DIASymGaussSeidelSmoother(
    const word &fieldName,
    const lduMatrix &matrix,
    const FieldField<Field, scalar> &interfaceBouCoeffs,
    const FieldField<Field, scalar> &interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList &interfaces)
    : lduMatrix::smoother(
          fieldName,
          matrix,
          interfaceBouCoeffs,
          interfaceIntCoeffs,
          interfaces),
      structuredMesh_(false),
      Nx_(0),
      Ny_(0),
      Nz_(0),
      useDIA_(false),
      upperCoeffI_(0.0),
      upperCoeffJ_(0.0),
      upperCoeffK_(0.0),
      lowerCoeffI_(0.0),
      lowerCoeffJ_(0.0),
      lowerCoeffK_(0.0),
      upperCoeffFieldI_(),
      upperCoeffFieldJ_(),
      upperCoeffFieldK_()
{
    const label nCells = matrix_.diag().size();

    const labelUList upperAddr = matrix_.lduAddr().upperAddr();
    const labelUList lowerAddr = matrix_.lduAddr().lowerAddr();
    label nFaces = upperAddr.size();
    const scalarField &upper = matrix_.upper();

    // unique offsets calculation

    std::set<label> offsets;
    for (label f = 0; f < lowerAddr.size(); f++)
    {
        offsets.insert(upperAddr[f] - lowerAddr[f]);
    }
    std::vector<label> sortedOffsets(offsets.begin(), offsets.end());

    // checking whether the size of offsets is 2 (for 2D) or 3 (for 3D)
    if (sortedOffsets.size() == 2)
    {
        if (sortedOffsets[0] == 1 && nCells % sortedOffsets[1] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = nCells / Nx_;
            Nz_ = 1;
            structuredMesh_ = true;
        }
    }
    else if (sortedOffsets.size() == 3)
    {
        if (sortedOffsets[0] == 1 && sortedOffsets[2] % sortedOffsets[1] == 0 && nCells % sortedOffsets[2] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = sortedOffsets[2] / sortedOffsets[1];
            Nz_ = nCells / sortedOffsets[2];
            structuredMesh_ = true;
        }
    }

    // // logging the outcome
    // // DEBUG: dumping first 5 faces per direction to see actual coefficient values

    // if(structuredMesh_)
    // {
    //     label countI = 0, countJ = 0, countK = 0;

    // }

    // checking coefficient uniformity
    bool allUniform = false;
    if (structuredMesh_ && !matrix_.asymmetric())
    {
        bool firstI = true, firstJ = true, firstK = true;

        allUniform = true;

        for (label f = 0; f < nFaces; f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            const scalar coeff = upper[f];

            if (offset == 1) // I
            {
                if (firstI)
                {
                    upperCoeffI_ = coeff;
                    firstI = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffI_), SMALL);
                    if (std::abs(coeff - upperCoeffI_) > tol)
                    {
                        allUniform = false;
                        break;
                    }
                }
            }
            else if (offset == Nx_) // J
            {
                if (firstJ)
                {
                    upperCoeffJ_ = coeff;
                    firstJ = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffJ_), SMALL);
                    if (std::abs(coeff - upperCoeffJ_) > tol)
                    {

                        allUniform = false;
                        break;
                    }
                }
            }
            else if (offset == Nx_ * Ny_) // K
            {
                if (firstK)
                {
                    upperCoeffK_ = coeff;
                    firstK = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffK_), SMALL);
                    if (std::abs(coeff - upperCoeffK_) > tol)
                    {

                        allUniform = false;
                        break;
                    }
                }
            }
        }
        if (allUniform)
        {
            lowerCoeffI_ = upperCoeffI_;
            lowerCoeffJ_ = upperCoeffJ_;
            lowerCoeffK_ = upperCoeffK_;
        }
    }

    useDIA_ = structuredMesh_ && !matrix.asymmetric();
    // useDIA_ = structuredMesh_ && !matrix.asymmetric() && (Nz_ > 1);

    if (useDIA_)
    {
        if (debug)
        {
            if (allUniform)
            {
                Info << "DIASymGaussSeidel: fast path enabled (uniform coefficients),  Nx=" << Nx_
                     << " Ny=" << Ny_ << " Nz=" << Nz_
                     << " upperCoeffs=(" << upperCoeffI_
                     << ", " << upperCoeffJ_
                     << ", " << upperCoeffK_ << ")"
                     << " (nCells=" << nCells << ")" << endl;
            }
            else
            {
                Info << "DIASymGaussSeidel: fast path enabled (variable coefficients),  Nx=" << Nx_
                     << " Ny=" << Ny_ << " Nz=" << Nz_
                     << " upperCoeffs=(" << upperCoeffI_
                     << ", " << upperCoeffJ_
                     << ", " << upperCoeffK_ << ")"
                     << " (nCells=" << nCells << ")" << endl;
            }
        }

        // initialising coefficient field
        upperCoeffFieldI_.setSize(nCells);
        upperCoeffFieldJ_.setSize(nCells);
        upperCoeffFieldK_.setSize(nCells);

        upperCoeffFieldI_ = 0.0;
        upperCoeffFieldJ_ = 0.0;
        upperCoeffFieldK_ = 0.0;

        for (label f = 0; f < nFaces; f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            const label owner = lowerAddr[f];

            if (offset == 1) // I
            {
                upperCoeffFieldI_[owner] = upper[f];
            }
            else if (offset == Nx_) // J
            {
                upperCoeffFieldJ_[owner] = upper[f];
            }
            else if (offset == Nx_ * Ny_) // K
            {
                upperCoeffFieldK_[owner] = upper[f];
            }
            else
            {
            }
        }
    }
    else if (!structuredMesh_)
    {
        if (debug)
        {
            Info << "DIASymGaussSeidel: mesh not structured, falling back (nCells="
                 << nCells << ")" << endl;
        }
    }
    else if (matrix_.asymmetric())
    {
        if (debug)
        {
            Info << "DIASymGaussSeidel: asymmetric matrix, falling back" << endl;
        }
    }

    if (useDIA_ && nCells == 125000) // only log on finest level of your test case
    {
        if (debug)
        {
            Info << "DEBUG: first 5 entries of upperCoeffFieldI_: ";
            for (label i = 0; i < 5; i++)
                Info << upperCoeffFieldI_[i] << " ";
            Info << endl;
            // same for J and K
            Info << "DEBUG: first 5 entries of upperCoeffFieldJ_: ";
            for (label i = 0; i < 5; i++)
                Info << upperCoeffFieldJ_[i] << " ";
            Info << endl;
            Info << "DEBUG: first 5 entries of upperCoeffFieldK_: ";
            for (label i = 0; i < 5; i++)
                Info << upperCoeffFieldK_[i] << " ";
            Info << endl;

            Info << "DEBUG: boundary slots (should be 0): "
                 << "upperCoeffFieldI_[" << (Nx_ - 1) << "]=" << upperCoeffFieldI_[Nx_ - 1] << " "
                 << "upperCoeffFieldJ_[" << (Nx_ * (Ny_ - 1)) << "]=" << upperCoeffFieldJ_[Nx_ * (Ny_ - 1)] << " "
                 << "upperCoeffFieldK_[" << (nCells - 1) << "]=" << upperCoeffFieldK_[nCells - 1]
                 << endl;
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DIASymGaussSeidelSmoother::smooth(
    scalarField &psi,
    const scalarField &source,
    const direction cmpt,
    const label nSweeps) const
{
    if (useDIA_)
    {
        scalar *__restrict__ psiPtr = psi.begin();
        const scalar *__restrict__ diagPtr = matrix_.diag().begin();
        const label nCells = psi.size();

        scalarField bPrime(nCells);
        scalar *__restrict__ bPrimePtr = bPrime.begin();

        // Cache Nx, Ny, Nz into local labels (avoids repeated member access in hot loop)
        label Nx = Nx_;
        label Ny = Ny_;
        label Nz = Nz_;

        const scalar *__restrict__ upperIPtr = upperCoeffFieldI_.begin();
        const scalar *__restrict__ upperJPtr = upperCoeffFieldJ_.begin();
        const scalar *__restrict__ upperKPtr = upperCoeffFieldK_.begin();

        const scalar *__restrict__ lowerIPtr = upperIPtr;
        const scalar *__restrict__ lowerJPtr = upperJPtr;
        const scalar *__restrict__ lowerKPtr = upperKPtr;

        label jStride = Nx;
        label kStride = Nx * Ny;

        // Parallel boundary initialisation.  The parallel boundary is treated
        // as an effective jacobi interface in the boundary.
        // Note: there is a change of sign in the coupled
        // interface update.  The reason for this is that the
        // internal coefficients are all located at the l.h.s. of
        // the matrix whereas the "implicit" coefficients on the
        // coupled boundaries are all created as if the
        // coefficient contribution is of a source-kind (i.e. they
        // have a sign as if they are on the r.h.s. of the matrix.
        // To compensate for this, it is necessary to turn the
        // sign of the contribution.

        FieldField<Field, scalar> &mBouCoeffs =
            const_cast<FieldField<Field, scalar> &>(
                interfaceBouCoeffs_);

        forAll(mBouCoeffs, patchi)
        {
            if (interfaces_.set(patchi))
            {
                mBouCoeffs[patchi].negate();
            }
        }

        // ---- Sweep Loop ----
        static double totalSweepTime = 0.0;
        static label totalCalls = 0;
        auto t0 = std::chrono::high_resolution_clock::now();
        LIKWID_MARKER_START("DIASYMsweep");
        for (label sweep = 0; sweep < nSweeps; sweep++)
        {
            bPrime = source;

            matrix_.initMatrixInterfaces(
                mBouCoeffs,
                interfaces_,
                psi,
                bPrime,
                cmpt);

            matrix_.updateMatrixInterfaces(
                mBouCoeffs,
                interfaces_,
                psi,
                bPrime,
                cmpt);

            // The DIA sweep
            for (label idx = 0; idx < nCells; idx++)
            {
                label j = idx % Nx;
                label i = (idx / Nx) % Ny;
                label k = idx / kStride;

                scalar psii = bPrimePtr[idx];

                // Upper triangle pull (reading forward neighbours)
                // Only happens if the forward neighbour exists (Not at boundary)
                if (j < Nx - 1)
                {
                    psii -= upperIPtr[idx] * psiPtr[idx + 1];
                }
                if (i < Ny - 1)
                {
                    psii -= upperJPtr[idx] * psiPtr[idx + jStride];
                }
                if (k < Nz - 1)
                {
                    psii -= upperKPtr[idx] * psiPtr[idx + kStride];
                }

                // divide by the diagnol
                psii /= diagPtr[idx];

                // Lower triangle push (updating forward neighbours' bPrime to be used in the next iteration)
                // Only happens if the forward neighbour exists (not at boundary)
                if (j < Nx - 1)
                {
                    bPrimePtr[idx + 1] -= lowerIPtr[idx] * psii;
                }
                if (i < Ny - 1)
                {
                    bPrimePtr[idx + jStride] -= lowerJPtr[idx] * psii;
                }
                if (k < Nz - 1)
                {
                    bPrimePtr[idx + kStride] -= lowerKPtr[idx] * psii;
                }

                psiPtr[idx] = psii;
            }

            for (label idx = nCells - 1; idx >= 0; idx--)
            {
                label j = idx % Nx;
                label i = (idx / Nx) % Ny;
                label k = idx / kStride;

                scalar psii = bPrimePtr[idx];

                if (j < Nx - 1)
                    psii -= upperIPtr[idx] * psiPtr[idx + 1];
                if (i < Ny - 1)
                    psii -= upperJPtr[idx] * psiPtr[idx + jStride];
                if (k < Nz - 1)
                    psii -= upperKPtr[idx] * psiPtr[idx + kStride];

                psii /= diagPtr[idx];

                if (j < Nx - 1)
                    bPrimePtr[idx + 1] -= lowerIPtr[idx] * psii;
                if (i < Ny - 1)
                    bPrimePtr[idx + jStride] -= lowerJPtr[idx] * psii;
                if (k < Nz - 1)
                    bPrimePtr[idx + kStride] -= lowerKPtr[idx] * psii;

                psiPtr[idx] = psii;
            }
        }
        LIKWID_MARKER_STOP("DIASYMsweep");
        auto t1 = std::chrono::high_resolution_clock::now();
        totalSweepTime += std::chrono::duration<double>(t1 - t0).count();
        totalCalls++;
        if (totalCalls % 1000 == 0)
        {
            Info << "DIASYMsweep: totalTime=" << totalSweepTime
                 << "s calls=" << totalCalls << endl;
        }

        // Restore interfaceBouCoeffs_
        forAll(mBouCoeffs, patchi)
        {
            if (interfaces_.set(patchi))
            {
                mBouCoeffs[patchi].negate();
            }
        }

        // LIKWID_MARKER_CLOSE;
    }
    else
    {
        Foam::symGaussSeidelSmoother::smooth(
            fieldName_,
            psi,
            matrix_,
            source,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt,
            nSweeps);
    }
}

// void Foam::DIASymGaussSeidelSmoother::smooth(
//     scalarField &psi,
//     const scalarField &source,
//     const direction cmpt,
//     const label nSweeps) const
// {
//     smooth(
//         fieldName_,
//         psi,
//         matrix_,
//         source,
//         interfaceBouCoeffs_,
//         interfaces_,
//         cmpt,
//         nSweeps);
// }

// ************************************************************************* //
