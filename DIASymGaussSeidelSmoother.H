/*---------------------------------------------------------------------------*\
Class
    Foam::DIAGaussSeidelSmoother

Description
    A lduMatrix::smoother for Gauss-Seidel, but works using the DIA
    addressing scheme, instead of the default LDU addressing scheme

SourceFiles
    DIAGaussSeidelSmoother.C
\*---------------------------------------------------------------------------*/

#ifndef DIAGaussSeidelSmoother_H
#define DIAGaussSeidelSmoother_H

#include "doubleScalar.H"
#include "lduMatrix.H"
#include "scalarField.H"

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                         Class DIAGaussSeidelSmoother Declaration
    \*---------------------------------------------------------------------------*/

    class DIAGaussSeidelSmoother : public lduMatrix::smoother
    {
    private:
        bool structuredMesh_;
        label Nx_;
        label Ny_;
        label Nz_;
        bool useDIA_;
        scalar upperCoeffI_;
        scalar upperCoeffJ_;
        scalar upperCoeffK_;
        scalar lowerCoeffI_;
        scalar lowerCoeffJ_;
        scalar lowerCoeffK_;
        scalarField upperCoeffFieldI_;
        scalarField upperCoeffFieldJ_;
        scalarField upperCoeffFieldK_;

    public:
        TypeName("DIAGaussSeidel");

        // Constructors

        //- Construct from components
        DIAGaussSeidelSmoother(
            const word &fieldName,
            const lduMatrix &matrix,
            const FieldField<Field, scalar> &interfaceBouCoeffs,
            const FieldField<Field, scalar> &interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList &interfaces);

        // Member Functions

        // // - Smooth for the given number of sweeps
        // static void smooth
        // (
        //     const word& fieldName,
        //     scalarField& psi,
        //     const lduMatrix& matrix,
        //     const scalarField& source,
        //     const FieldField<Field, scalar>& interfaceBouCoeffs,
        //     const lduInterfaceFieldPtrsList& interfaces,
        //     const direction cmpt,
        //     const label nSweeps
        // );

        //- Smooth the solution for a given number of sweeps
        virtual void smooth(
            scalarField &psi,
            const scalarField &Source,
            const direction cmpt,
            const label nSweeps) const;
    };
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
