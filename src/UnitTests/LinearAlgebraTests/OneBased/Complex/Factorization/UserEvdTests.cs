﻿// <copyright file="UserEvdTests.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
// Copyright (c) 2009-2015 Math.NET
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using MathNet.Numerics.LinearAlgebra.OneBased;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Complex.Factorization
{
#if NOSYSNUMERICS
    using Complex = Numerics.Complex;
#else
    using Complex = System.Numerics.Complex;
#endif

    /// <summary>
    /// Eigenvalues factorization tests for an user matrix.
    /// </summary>
    [TestFixture, Category("LAFactorization")]
    public class UserEvdTests
    {
        /// <summary>
        /// Can factorize identity matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(10)]
        [TestCase(100)]
        public void CanFactorizeIdentity(int order)
        {
            var matrixI = UserDefinedMatrix.Identity(order);
            var factorEvd = matrixI.Evd();
            var eigenValues = factorEvd.EigenValues;
            var eigenVectors = factorEvd.EigenVectors;
            var d = factorEvd.D;

            Assert.AreEqual(matrixI.RowCount, eigenVectors.RowCount);
            Assert.AreEqual(matrixI.RowCount, eigenVectors.ColumnCount);

            Assert.AreEqual(matrixI.ColumnCount, d.RowCount);
            Assert.AreEqual(matrixI.ColumnCount, d.ColumnCount);
            AssertHelpers.VectorHasValue(eigenValues, Complex.One);
        }

        /// <summary>
        /// Can factorize a random square matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanFactorizeRandomMatrix(int order)
        {
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.Random(order, order, 1).ToArray());
            var factorEvd = matrixA.Evd();
            var eigenVectors = factorEvd.EigenVectors;
            var d = factorEvd.D;

            Assert.AreEqual(order, eigenVectors.RowCount);
            Assert.AreEqual(order, eigenVectors.ColumnCount);

            Assert.AreEqual(order, d.RowCount);
            Assert.AreEqual(order, d.ColumnCount);

            // Make sure the A*V = λ*V 
            var matrixAv = matrixA * eigenVectors;
            var matrixLv = eigenVectors * d;
            AssertHelpers.AlmostEqualRelative(matrixAv, matrixLv, 7);
        }

        /// <summary>
        /// Can factorize a symmetric random square matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanFactorizeRandomSymmetricMatrix(int order)
        {
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            var factorEvd = matrixA.Evd(Symmetricity.Hermitian);
            var eigenVectors = factorEvd.EigenVectors;
            var d = factorEvd.D;

            Assert.AreEqual(order, eigenVectors.RowCount);
            Assert.AreEqual(order, eigenVectors.ColumnCount);

            Assert.AreEqual(order, d.RowCount);
            Assert.AreEqual(order, d.ColumnCount);

            // Make sure the A = V*λ*VT 
            var matrix = eigenVectors * d * eigenVectors.ConjugateTranspose();
            AssertHelpers.AlmostEqualRelative(matrix, matrixA, 7);
        }

        /// <summary>
        /// Can check rank of square matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanCheckRankSquare(int order)
        {
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.Random(order, order, 1).ToArray());
            var factorEvd = matrixA.Evd();

            Assert.AreEqual(factorEvd.Rank, order);
        }

        /// <summary>
        /// Can check rank of square singular matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanCheckRankOfSquareSingular(int order)
        {
            var matrixA = new UserDefinedMatrix(order, order);
            matrixA[1, 1] = 1;
            matrixA[order, order] = 1;
            for (var i = 2; i < order; i++)
            {
                matrixA[i, i - 1] = 1;
                matrixA[i, i + 1] = 1;
                matrixA[i - 1, i] = 1;
                matrixA[i + 1, i] = 1;
            }

            var factorEvd = matrixA.Evd();

            Assert.AreEqual(factorEvd.Determinant, Complex.Zero);
            Assert.AreEqual(factorEvd.Rank, order - 1);
        }

        /// <summary>
        /// Identity determinant is one.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(10)]
        [TestCase(100)]
        public void IdentityDeterminantIsOne(int order)
        {
            var matrixI = UserDefinedMatrix.Identity(order);
            var factorEvd = matrixI.Evd();
            Assert.AreEqual(Complex.One, factorEvd.Determinant);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random vector and symmetric matrix (Ax=b).
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [Test]
        public void CanSolveForRandomVectorAndSymmetricMatrix([Values(1, 2, 5, 10, 50, 100)] int order)
        {
            var A = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            MatrixHelpers.ForceHermitian(A);
            var ACopy = A.Clone();
            var evd = A.Evd(Symmetricity.Hermitian);

            var b = new UserDefinedVector(Vector<Complex>.Build.Random(order, 1).ToArray());
            var bCopy = b.Clone();

            var x = evd.Solve(b);

            var bReconstruct = A * x;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(b, bReconstruct, 9);

            // Make sure A/B didn't change.
            AssertHelpers.AlmostEqual(ACopy, A, 14);
            AssertHelpers.AlmostEqual(bCopy, b, 14);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix and symmetric matrix (AX=B).
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [Test]
        public void CanSolveForRandomMatrixAndSymmetricMatrix([Values(1, 2, 5, 10, 50, 100)] int order)
        {
            var A = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            MatrixHelpers.ForceHermitian(A);
            var ACopy = A.Clone();
            var evd = A.Evd(Symmetricity.Hermitian);

            var B = new UserDefinedMatrix(Matrix<Complex>.Build.Random(order, order, 1).ToArray());
            var BCopy = B.Clone();

            var X = evd.Solve(B);

            // The solution X row dimension is equal to the column dimension of A
            Assert.AreEqual(A.ColumnCount, X.RowCount);

            // The solution X has the same number of columns as B
            Assert.AreEqual(B.ColumnCount, X.ColumnCount);

            var BReconstruct = A * X;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(B, BReconstruct, 9);

            // Make sure A/B didn't change.
            AssertHelpers.AlmostEqual(ACopy, A, 14);
            AssertHelpers.AlmostEqual(BCopy, B, 14);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random vector and symmetric matrix (Ax=b) into a result matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [Test]
        public void CanSolveForRandomVectorAndSymmetricMatrixWhenResultVectorGiven([Values(1, 2, 5, 10, 50, 100)] int order)
        {
            var A = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            MatrixHelpers.ForceHermitian(A);
            var ACopy = A.Clone();
            var evd = A.Evd(Symmetricity.Hermitian);

            var b = new UserDefinedVector(Vector<Complex>.Build.Random(order, 1).ToArray());
            var bCopy = b.Clone();

            var x = new UserDefinedVector(order);
            evd.Solve(b, x);

            var bReconstruct = A * x;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(b, bReconstruct, 9);

            // Make sure A/B didn't change.
            AssertHelpers.AlmostEqual(ACopy, A, 14);
            AssertHelpers.AlmostEqual(bCopy, b, 14);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix and symmetric matrix (AX=B) into result matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [Test]
        public void CanSolveForRandomMatrixAndSymmetricMatrixWhenResultMatrixGiven([Values(1, 2, 5, 10, 50, 100)] int order)
        {
            var A = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            MatrixHelpers.ForceHermitian(A);
            var ACopy = A.Clone();
            var evd = A.Evd(Symmetricity.Hermitian);

            var B = new UserDefinedMatrix(Matrix<Complex>.Build.Random(order, order, 1).ToArray());
            var BCopy = B.Clone();

            var X = new UserDefinedMatrix(order, order);
            evd.Solve(B, X);

            // The solution X row dimension is equal to the column dimension of A
            Assert.AreEqual(A.ColumnCount, X.RowCount);

            // The solution X has the same number of columns as B
            Assert.AreEqual(B.ColumnCount, X.ColumnCount);

            var BReconstruct = A * X;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(B, BReconstruct, 9);

            // Make sure A/B didn't change.
            AssertHelpers.AlmostEqual(ACopy, A, 14);
            AssertHelpers.AlmostEqual(BCopy, B, 14);
        }
    }
}