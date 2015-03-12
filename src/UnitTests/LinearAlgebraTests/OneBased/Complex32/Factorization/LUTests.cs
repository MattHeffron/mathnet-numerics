// <copyright file="LUTests.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.OneBased.Complex32;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Complex32.Factorization
{
    using Numerics;

    /// <summary>
    /// LU factorization tests for a dense matrix.
    /// </summary>
    [TestFixture, Category("LAFactorization")]
    public class LUTests
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
            var matrixI = DenseMatrix.CreateIdentity(order);
            var factorLU = matrixI.LU();

            // Check lower triangular part.
            var matrixL = factorLU.L;
            Assert.AreEqual(matrixI.RowCount, matrixL.RowCount);
            Assert.AreEqual(matrixI.ColumnCount, matrixL.ColumnCount);
            AssertHelpers.IsDiagonal(matrixL);
            AssertHelpers.DiagonalHasValue(matrixL, Complex32.One);

            // Check upper triangular part.
            var matrixU = factorLU.U;
            Assert.AreEqual(matrixI.RowCount, matrixU.RowCount);
            Assert.AreEqual(matrixI.ColumnCount, matrixU.ColumnCount);
            AssertHelpers.IsDiagonal(matrixU);
            AssertHelpers.DiagonalHasValue(matrixU, Complex32.One);
        }

        /// <summary>
        /// LU factorization fails with a non-square matrix.
        /// </summary>
        [Test]
        public void LUFailsWithNonSquareMatrix()
        {
            var matrix = new DenseMatrix(3, 1);
            Assert.That(() => matrix.LU(), Throws.ArgumentException);
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
            var matrixI = DenseMatrix.CreateIdentity(order);
            var lu = matrixI.LU();
            Assert.AreEqual(Complex32.One, lu.Determinant);
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
            var matrixX = Matrix<Complex32>.Build.Random(order, order, 1);
            var factorLU = matrixX.LU();
            var matrixL = factorLU.L;
            var matrixU = factorLU.U;

            // Make sure the factors have the right dimensions.
            Assert.AreEqual(order, matrixL.RowCount);
            Assert.AreEqual(order, matrixL.ColumnCount);
            Assert.AreEqual(order, matrixU.RowCount);
            Assert.AreEqual(order, matrixU.ColumnCount);

            // Make sure the L factor is lower triangular.
            AssertHelpers.IsLowerTriangular(matrixL);

            // Make sure the U factor is upper triangular.
            AssertHelpers.IsUpperTriangular(matrixU);

            // Make sure the LU factor times it's transpose is the original matrix.
            var matrixXfromLU = matrixL * matrixU;
            var permutationInverse = factorLU.P.Inverse();
            matrixXfromLU.PermuteRows(permutationInverse);
            AssertHelpers.AlmostEqualRelative(matrixX, matrixXfromLU, 3);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random vector (Ax=b).
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanSolveForRandomVector(int order)
        {
            var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixACopy = matrixA.Clone();
            var factorLU = matrixA.LU();

            var vectorb = Vector<Complex32>.Build.Random(order, 1);
            var resultx = factorLU.Solve(vectorb);

            Assert.AreEqual(matrixA.ColumnCount, resultx.Count);

            var matrixBReconstruct = matrixA * resultx;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(vectorb, matrixBReconstruct, 3);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix (AX=B).
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanSolveForRandomMatrix(int order)
        {
            var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixACopy = matrixA.Clone();
            var factorLU = matrixA.LU();

            var matrixB = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixX = factorLU.Solve(matrixB);

            // The solution X row dimension is equal to the column dimension of A
            Assert.AreEqual(matrixA.ColumnCount, matrixX.RowCount);

            // The solution X has the same number of columns as B
            Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

            var matrixBReconstruct = matrixA * matrixX;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(matrixB, matrixBReconstruct, 3);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);
        }

        /// <summary>
        /// Can solve for a random vector into a result vector.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanSolveForRandomVectorWhenResultVectorGiven(int order)
        {
            var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixACopy = matrixA.Clone();
            var factorLU = matrixA.LU();
            var vectorb = Vector<Complex32>.Build.Random(order, 1);
            var vectorbCopy = vectorb.Clone();
            var resultx = new DenseVector(order);
            factorLU.Solve(vectorb, resultx);

            Assert.AreEqual(vectorb.Count, resultx.Count);

            var matrixBReconstruct = matrixA * resultx;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(vectorb, matrixBReconstruct, 3);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);

            // Make sure b didn't change.
            AssertHelpers.AreEqual(vectorbCopy, vectorb);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix (AX=B) into a result matrix.
        /// </summary>
        /// <param name="order">Matrix row number.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanSolveForRandomMatrixWhenResultMatrixGiven(int order)
        {
            var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixACopy = matrixA.Clone();
            var factorLU = matrixA.LU();

            var matrixB = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixBCopy = matrixB.Clone();

            var matrixX = new DenseMatrix(order, order);
            factorLU.Solve(matrixB, matrixX);

            // The solution X row dimension is equal to the column dimension of A
            Assert.AreEqual(matrixA.ColumnCount, matrixX.RowCount);

            // The solution X has the same number of columns as B
            Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

            var matrixBReconstruct = matrixA * matrixX;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(matrixB, matrixBReconstruct, 3);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);

            // Make sure B didn't change.
            AssertHelpers.AreEqual(matrixBCopy, matrixB);
        }

        /// <summary>
        /// Can inverse a matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(1)]
        [TestCase(2)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(50)]
        [TestCase(100)]
        public void CanInverse(int order)
        {
            var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixACopy = matrixA.Clone();
            var factorLU = matrixA.LU();

            var matrixAInverse = factorLU.Inverse();

            // The inverse dimension is equal A
            Assert.AreEqual(matrixAInverse.RowCount, matrixAInverse.RowCount);
            Assert.AreEqual(matrixAInverse.ColumnCount, matrixAInverse.ColumnCount);

            var matrixIdentity = matrixA * matrixAInverse;

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);

            // Check if multiplication of A and AI produced identity matrix.
            AssertHelpers.ValuesAssertion(matrixIdentity, (i, j, v) => AssertHelpers.AlmostEqualRelative(i == j ? Complex.One : Complex.Zero, matrixIdentity[i, j], 3));
        }
    }
}
