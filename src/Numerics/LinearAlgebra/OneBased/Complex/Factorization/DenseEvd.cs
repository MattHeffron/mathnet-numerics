// <copyright file="DenseEvd.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
//
// Copyright (c) 2009-2015 Math.NET
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;
using MathNet.Numerics.Properties;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Complex.Factorization
{

#if NOSYSNUMERICS
    using Numerics;
#else
    using System.Numerics;
#endif

    /// <summary>
    /// Eigenvalues and eigenvectors of a complex matrix.
    /// </summary>
    /// <remarks>
    /// If A is hermitan, then A = V*D*V' where the eigenvalue matrix D is
    /// diagonal and the eigenvector matrix V is hermitan.
    /// I.e. A = V*D*V' and V*VH=I.
    /// If A is not symmetric, then the eigenvalue matrix D is block diagonal
    /// with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    /// lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    /// columns of V represent the eigenvectors in the sense that A*V = V*D,
    /// i.e. A.Multiply(V) equals V.Multiply(D).  The matrix V may be badly
    /// conditioned, or even singular, so the validity of the equation
    /// A = V*D*Inverse(V) depends upon V.Condition().
    /// </remarks>
    internal sealed class DenseEvd : Evd
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="DenseEvd"/> class. This object will compute the
        /// the eigenvalue decomposition when the constructor is called and cache it's decomposition.
        /// </summary>
        /// <param name="matrix">The matrix to factor.</param>
        /// <param name="symmetricity">If it is known whether the matrix is symmetric or not the routine can skip checking it itself.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="matrix"/> is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">If EVD algorithm failed to converge with matrix <paramref name="matrix"/>.</exception>
        public static DenseEvd Create(DenseMatrix matrix, Symmetricity symmetricity)
        {
            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            var order = matrix.RowCount;

            // Initialize matrices for eigenvalues and eigenvectors
            var eigenVectors = DenseMatrix.CreateIdentity(order);
            var blockDiagonal = new DenseMatrix(order);
            var eigenValues = new DenseVector(order);

            bool isSymmetric;
            switch (symmetricity)
            {
                case Symmetricity.Hermitian:
                    isSymmetric = true;
                    break;
                case Symmetricity.Asymmetric:
                    isSymmetric = false;
                    break;
                default:
                    isSymmetric = matrix.IsHermitian();
                    break;
            }

            Control.LinearAlgebraProvider.EigenDecomp(isSymmetric, order, matrix.Values, eigenVectors.Values, eigenValues.Values, blockDiagonal.Values);

            return new DenseEvd(eigenVectors, eigenValues, blockDiagonal, isSymmetric);
        }

        DenseEvd(Matrix<Complex> eigenVectors, Vector<Complex> eigenValues, Matrix<Complex> blockDiagonal, bool isSymmetric)
            : base(eigenVectors, eigenValues, blockDiagonal, isSymmetric)
        {
        }

        /// <summary>
        /// Reduces a complex hermitian matrix to a real symmetric tridiagonal matrix using unitary similarity transformations.
        /// </summary>
        /// <param name="matrixA">Source matrix to reduce</param>
        /// <param name="d">Output: Arrays for internal storage of real parts of eigenvalues</param>
        /// <param name="e">Output: Arrays for internal storage of imaginary parts of eigenvalues</param>
        /// <param name="tau">Output: Arrays that contains further information about the transformations.</param>
        /// <param name="order">Order of initial matrix</param>
        /// <remarks>This is derived from the Algol procedures HTRIDI by 
        /// Smith, Boyle, Dongarra, Garbow, Ikebe, Klema, Moler, and Wilkinson, Handbook for 
        /// Auto. Comp., Vol.ii-Linear Algebra, and the corresponding 
        /// Fortran subroutine in EISPACK.</remarks>
        internal static void SymmetricTridiagonalize(Complex[] matrixA, double[] d, double[] e, Complex[] tau, int order)
        {
            // Since there is no dependence on the zero-based/one-based differences for this, use the same implementation for both
            // This protects against possible version skew if it is ever changed.
            LinearAlgebra.Complex.Factorization.DenseEvd.SymmetricTridiagonalize(matrixA, d, e, tau, order);
        }

        /// <summary>
        /// Symmetric tridiagonal QL algorithm.
        /// </summary>
        /// <param name="dataEv">Data array of matrix V (eigenvectors)</param>
        /// <param name="d">Arrays for internal storage of real parts of eigenvalues</param>
        /// <param name="e">Arrays for internal storage of imaginary parts of eigenvalues</param>
        /// <param name="order">Order of initial matrix</param>
        /// <remarks>This is derived from the Algol procedures tql2, by
        /// Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        /// Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        /// Fortran subroutine in EISPACK.</remarks>
        /// <exception cref="NonConvergenceException"></exception>
        internal static void SymmetricDiagonalize(Complex[] dataEv, double[] d, double[] e, int order)
        {
            // Since there is no dependence on the zero-based/one-based differences for this, use the same implementation for both
            // This protects against possible version skew if it is ever changed.
            LinearAlgebra.Complex.Factorization.DenseEvd.SymmetricDiagonalize(dataEv, d, e, order);
        }

        /// <summary>
        /// Determines eigenvectors by undoing the symmetric tridiagonalize transformation
        /// </summary>
        /// <param name="dataEv">Data array of matrix V (eigenvectors)</param>
        /// <param name="matrixA">Previously tridiagonalized matrix by <see cref="SymmetricTridiagonalize"/>.</param>
        /// <param name="tau">Contains further information about the transformations</param>
        /// <param name="order">Input matrix order</param>
        /// <remarks>This is derived from the Algol procedures HTRIBK, by
        /// by Smith, Boyle, Dongarra, Garbow, Ikebe, Klema, Moler, and Wilkinson, Handbook for
        /// Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        /// Fortran subroutine in EISPACK.</remarks>
        internal static void SymmetricUntridiagonalize(Complex[] dataEv, Complex[] matrixA, Complex[] tau, int order)
        {
            // Since there is no dependence on the zero-based/one-based differences for this, use the same implementation for both
            // This protects against possible version skew if it is ever changed.
            LinearAlgebra.Complex.Factorization.DenseEvd.SymmetricUntridiagonalize(dataEv, matrixA, tau, order);
        }

        /// <summary>
        /// Nonsymmetric reduction to Hessenberg form.
        /// </summary>
        /// <param name="dataEv">Data array of matrix V (eigenvectors)</param>
        /// <param name="matrixH">Array for internal storage of nonsymmetric Hessenberg form.</param>
        /// <param name="order">Order of initial matrix</param>
        /// <remarks>This is derived from the Algol procedures orthes and ortran,
        /// by Martin and Wilkinson, Handbook for Auto. Comp.,
        /// Vol.ii-Linear Algebra, and the corresponding
        /// Fortran subroutines in EISPACK.</remarks>
        internal static void NonsymmetricReduceToHessenberg(Complex[] dataEv, Complex[] matrixH, int order)
        {
            // Since there is no dependence on the zero-based/one-based differences for this, use the same implementation for both
            // This protects against possible version skew if it is ever changed.
            LinearAlgebra.Complex.Factorization.DenseEvd.NonsymmetricReduceToHessenberg(dataEv, matrixH, order);
        }

        /// <summary>
        /// Nonsymmetric reduction from Hessenberg to real Schur form.
        /// </summary>
        /// <param name="vectorV">Data array of the eigenvectors</param>
        /// <param name="dataEv">Data array of matrix V (eigenvectors)</param>
        /// <param name="matrixH">Array for internal storage of nonsymmetric Hessenberg form.</param>
        /// <param name="order">Order of initial matrix</param>
        /// <remarks>This is derived from the Algol procedure hqr2,
        /// by Martin and Wilkinson, Handbook for Auto. Comp.,
        /// Vol.ii-Linear Algebra, and the corresponding
        /// Fortran subroutine in EISPACK.</remarks>
        internal static void NonsymmetricReduceHessenberToRealSchur(Complex[] vectorV, Complex[] dataEv, Complex[] matrixH, int order)
        {
            // Since there is no dependence on the zero-based/one-based differences for this, use the same implementation for both
            // This protects against possible version skew if it is ever changed.
            LinearAlgebra.Complex.Factorization.DenseEvd.NonsymmetricReduceHessenberToRealSchur(vectorV, dataEv, matrixH, order);
        }

        /// <summary>
        /// Solves a system of linear equations, <b>AX = B</b>, with A SVD factorized.
        /// </summary>
        /// <param name="input">The right hand side <see cref="Matrix{T}"/>, <b>B</b>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <b>X</b>.</param>
        public override void Solve(Matrix<Complex> input, Matrix<Complex> result)
        {
            // The solution X should have the same number of columns as B
            if (input.ColumnCount != result.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            }

            // The dimension compatibility conditions for X = A\B require the two matrices A and B to have the same number of rows
            if (EigenValues.Count != input.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameRowDimension);
            }

            // The solution X row dimension is equal to the column dimension of A
            if (EigenValues.Count != result.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            }

            if (IsSymmetric)
            {
                var order = EigenValues.Count;
                var tmp = new Complex[order];

                for (var k = 1; k <= order; k++)
                {
                    for (int j = 1, j0 = 0; j <= order; j++, j0++)
                    {
                        Complex value = 0.0;
                        if (j <= order)     // CONSIDER: Isn't j *always* <= order by construction??
                        {
                            for (var i = 0; i < order; i++)
                            {
                                value += ((DenseMatrix) EigenVectors).Values[(j0*order) + i].Conjugate()*input.At(i + 1, k);
                            }

                            value /= EigenValues[j].Real;
                        }

                        tmp[j0] = value;
                    }

                    for (int j = 1, j0 = 0; j <= order; j++, j0++)
                    {
                        Complex value = 0.0;
                        for (var i = 0; i < order; i++)
                        {
                            value += ((DenseMatrix) EigenVectors).Values[(i*order) + j0]*tmp[i];
                        }

                        result.At(j, k, value);
                    }
                }
            }
            else
            {
                throw new ArgumentException(Resources.ArgumentMatrixSymmetric);
            }
        }

        /// <summary>
        /// Solves a system of linear equations, <b>Ax = b</b>, with A EVD factorized.
        /// </summary>
        /// <param name="input">The right hand side vector, <b>b</b>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <b>x</b>.</param>
        public override void Solve(Vector<Complex> input, Vector<Complex> result)
        {
            // Ax=b where A is an m x m matrix
            // Check that b is a column vector with m entries
            if (EigenValues.Count != input.Count)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            // Check that x is a column vector with n entries
            if (EigenValues.Count != result.Count)
            {
                throw Matrix.DimensionsDontMatch<ArgumentException>(EigenValues, result);
            }

            if (IsSymmetric)
            {
                // Symmetric case -> x = V * inv(λ) * VH * b;
                var order = EigenValues.Count;
                var tmp = new Complex[order];
                Complex value;

                for (int j = 1, j0 = 0; j <= order; j++, j0++)
                {
                    value = 0;
                    if (j <= order)     // CONSIDER: Isn't j *always* <= order by construction??
                    {
                        for (var i = 0; i < order; i++)
                        {
                            value += ((DenseMatrix) EigenVectors).Values[(j0*order) + i].Conjugate()*input.At(i + 1);
                        }

                        value /= EigenValues[j].Real;
                    }

                    tmp[j0] = value;
                }

                for (int j = 1, j0 = 0; j <= order; j++, j0++)
                {
                    value = 0;
                    for (var i = 0; i < order; i++)
                    {
                        value += ((DenseMatrix) EigenVectors).Values[(i*order) + j0]*tmp[i];
                    }

                    result.At(j, value);
                }
            }
            else
            {
                throw new ArgumentException(Resources.ArgumentMatrixSymmetric);
            }
        }
    }
}
