// <copyright file="UserCholeskyTests.cs" company="Math.NET">
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
    /// Cholesky factorization tests for a user matrix.
    /// </summary>
    [TestFixture, Category("LAFactorization")]
    public class UserCholeskyTests
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
            var factorC = matrixI.Cholesky().Factor;

            Assert.AreEqual(matrixI.RowCount, factorC.RowCount);
            Assert.AreEqual(matrixI.ColumnCount, factorC.ColumnCount);

            AssertHelpers.IsIdentity(factorC);
        }

        /// <summary>
        /// Cholesky factorization fails with diagonal a non-positive definite matrix.
        /// </summary>
        [Test]
        public void CholeskyFailsWithDiagonalNonPositiveDefiniteMatrix()
        {
            var matrixI = UserDefinedMatrix.Identity(8);
            matrixI[4, 4] = -4.0;
            Assert.That(() => matrixI.Cholesky(), Throws.ArgumentException);
        }

        /// <summary>
        /// Cholesky factorization fails with a non-square matrix.
        /// </summary>
        [Test]
        public void CholeskyFailsWithNonSquareMatrix()
        {
            var matrixI = new UserDefinedMatrix(3, 2);
            Assert.That(() => matrixI.Cholesky(), Throws.ArgumentException);
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
            var factorC = matrixI.Cholesky();
            Assert.AreEqual(Complex.One, factorC.Determinant);
            Assert.AreEqual(Complex.Zero, factorC.DeterminantLn);
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
            var matrixX = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            var chol = matrixX.Cholesky();
            var factorC = chol.Factor;

            // Make sure the Cholesky factor has the right dimensions.
            Assert.AreEqual(order, factorC.RowCount);
            Assert.AreEqual(order, factorC.ColumnCount);

            // Make sure the Cholesky factor is lower triangular.
            AssertHelpers.IsLowerTriangular(factorC);

            // Make sure the cholesky factor times it's transpose is the original matrix.
            var matrixXfromC = factorC * factorC.ConjugateTranspose();
            AssertHelpers.AlmostEqualRelative(matrixX, matrixXfromC, 8);
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
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            var matrixACopy = matrixA.Clone();
            var chol = matrixA.Cholesky();
            var b = new UserDefinedVector(Vector<Complex>.Build.Random(order, 1).ToArray());
            var x = chol.Solve(b);

            Assert.AreEqual(b.Count, x.Count);

            var matrixBReconstruct = matrixA * x;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(b, matrixBReconstruct, 10);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix (AX=B).
        /// </summary>
        /// <param name="row">Matrix row number.</param>
        /// <param name="col">Matrix column number.</param>
        [TestCase(1, 1)]
        [TestCase(2, 4)]
        [TestCase(5, 8)]
        [TestCase(10, 3)]
        [TestCase(50, 10)]
        [TestCase(100, 100)]
        public void CanSolveForRandomMatrix(int row, int col)
        {
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(row, 1).ToArray());
            var matrixACopy = matrixA.Clone();
            var chol = matrixA.Cholesky();
            var matrixB = new UserDefinedMatrix(Matrix<Complex>.Build.Random(row, col, 1).ToArray());
            var matrixX = chol.Solve(matrixB);

            Assert.AreEqual(matrixB.RowCount, matrixX.RowCount);
            Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

            var matrixBReconstruct = matrixA * matrixX;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(matrixB, matrixBReconstruct, 10);

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
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(order, 1).ToArray());
            var matrixACopy = matrixA.Clone();
            var chol = matrixA.Cholesky();
            var b = new UserDefinedVector(Vector<Complex>.Build.Random(order, 1).ToArray());
            var matrixBCopy = b.Clone();
            var x = new UserDefinedVector(order);
            chol.Solve(b, x);

            Assert.AreEqual(b.Count, x.Count);

            var matrixBReconstruct = matrixA * x;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(b, matrixBReconstruct, 10);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);

            // Make sure b didn't change.
            AssertHelpers.AreEqual(matrixBCopy, b);
        }

        /// <summary>
        /// Can solve a system of linear equations for a random matrix (AX=B) into a result matrix.
        /// </summary>
        /// <param name="row">Matrix row number.</param>
        /// <param name="col">Matrix column number.</param>
        [TestCase(1, 1)]
        [TestCase(2, 4)]
        [TestCase(5, 8)]
        [TestCase(10, 3)]
        [TestCase(50, 10)]
        [TestCase(100, 100)]
        public void CanSolveForRandomMatrixWhenResultMatrixGiven(int row, int col)
        {
            var matrixA = new UserDefinedMatrix(Matrix<Complex>.Build.RandomPositiveDefinite(row, 1).ToArray());
            var matrixACopy = matrixA.Clone();
            var chol = matrixA.Cholesky();
            var matrixB = new UserDefinedMatrix(Matrix<Complex>.Build.Random(row, col, 1).ToArray());
            var matrixBCopy = matrixB.Clone();
            var matrixX = new UserDefinedMatrix(row, col);
            chol.Solve(matrixB, matrixX);

            Assert.AreEqual(matrixB.RowCount, matrixX.RowCount);
            Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

            var matrixBReconstruct = matrixA * matrixX;

            // Check the reconstruction.
            AssertHelpers.AlmostEqual(matrixB, matrixBReconstruct, 10);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);

            // Make sure B didn't change.
            AssertHelpers.AreEqual(matrixBCopy, matrixB);
        }
    }
}
