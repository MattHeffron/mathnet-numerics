// <copyright file="GramSchmidtTests.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.OneBased.Complex32.Factorization;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Complex32.Factorization
{
    using Numerics;

    /// <summary>
    /// GramSchmidt factorization tests for a dense matrix.
    /// </summary>
    [TestFixture, Category("LAFactorization")]
    public class GramSchmidtTests
    {
        /// <summary>
        /// Constructor with wide matrix throws <c>ArgumentException</c>.
        /// </summary>
        [Test]
        public void ConstructorWideMatrixThrowsInvalidMatrixOperationException()
        {
            Assert.That(() => DenseGramSchmidt.Create(new DenseMatrix(3, 4)), Throws.ArgumentException);
        }

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
            var factorGramSchmidt = matrixI.GramSchmidt();
            var q = factorGramSchmidt.Q;
            var r = factorGramSchmidt.R;

            Assert.AreEqual(matrixI.RowCount, q.RowCount);
            Assert.AreEqual(matrixI.ColumnCount, q.ColumnCount);
            AssertHelpers.IsDiagonal(r);
            AssertHelpers.DiagonalHasValue(r, Complex32.One);
            AssertHelpers.IsDiagonal(q);
            AssertHelpers.DiagonalHasValue(q, Complex32.One);
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
            var factorGramSchmidt = matrixI.GramSchmidt();
            Assert.AreEqual(Complex32.One, factorGramSchmidt.Determinant);
        }

        /// <summary>
        /// Can factorize a random matrix.
        /// </summary>
        /// <param name="row">Matrix row number.</param>
        /// <param name="column">Matrix column number.</param>
        [TestCase(1, 1)]
        [TestCase(2, 2)]
        [TestCase(5, 5)]
        [TestCase(10, 6)]
        [TestCase(50, 48)]
        [TestCase(100, 98)]
        public void CanFactorizeRandomMatrix(int row, int column)
        {
            var matrixA = Matrix<Complex32>.Build.Random(row, column, 1);
            var factorGramSchmidt = matrixA.GramSchmidt();
            var q = factorGramSchmidt.Q;
            var r = factorGramSchmidt.R;

            // Make sure the Q has the right dimensions.
            Assert.AreEqual(row, q.RowCount);
            Assert.AreEqual(column, q.ColumnCount);

            // Make sure the R has the right dimensions.
            Assert.AreEqual(column, r.RowCount);
            Assert.AreEqual(column, r.ColumnCount);

            // Make sure the R factor is upper triangular.
            AssertHelpers.IsUpperTriangular(r);

            // Make sure the Q*R is the original matrix.
            var matrixQfromR = q * r;
            AssertHelpers.AlmostEqualRelative(matrixA, matrixQfromR, 3);

            // Make sure the Q is unitary --> (Q*)x(Q) = I
            var matrixQсtQ = q.ConjugateTranspose() * q;
            AssertHelpers.IsDiagonal(matrixQсtQ);
            AssertHelpers.DiagonalHasValue(matrixQсtQ, Complex32.One);
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
            var factorGramSchmidt = matrixA.GramSchmidt();

            var vectorb = Vector<Complex32>.Build.Random(order, 1);
            var resultx = factorGramSchmidt.Solve(vectorb);

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
            var factorGramSchmidt = matrixA.GramSchmidt();

            var matrixB = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixX = factorGramSchmidt.Solve(matrixB);

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
            var factorGramSchmidt = matrixA.GramSchmidt();
            var vectorb = Vector<Complex32>.Build.Random(order, 1);
            var vectorbCopy = vectorb.Clone();
            var resultx = new DenseVector(order);
            factorGramSchmidt.Solve(vectorb, resultx);

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
        /// <param name="order">Matrix order.</param>
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
            var factorGramSchmidt = matrixA.GramSchmidt();

            var matrixB = Matrix<Complex32>.Build.Random(order, order, 1);
            var matrixBCopy = matrixB.Clone();

            var matrixX = new DenseMatrix(order, order);
            factorGramSchmidt.Solve(matrixB, matrixX);

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
        /// Can solve when using a tall matrix.
        /// </summary>
        [Test]
        public void CanSolveForMatrixWithTallRandomMatrix()
        {
            var matrixA = Matrix<Complex32>.Build.Random(20, 10, 1);
            var matrixACopy = matrixA.Clone();
            var factorQR = matrixA.GramSchmidt();

            var matrixB = Matrix<Complex32>.Build.Random(20, 5, 1);
            var matrixX = factorQR.Solve(matrixB);

            // The solution X row dimension is equal to the column dimension of A
            Assert.AreEqual(matrixA.ColumnCount, matrixX.RowCount);

            // The solution X has the same number of columns as B
            Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

            var test = (matrixA.ConjugateTranspose() * matrixA).Inverse() * matrixA.ConjugateTranspose() * matrixB;

            AssertHelpers.AlmostEqual(test, matrixX, 5);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);
        }

        /// <summary>
        /// Can solve when using a tall matrix.
        /// </summary>
        [Test]
        public void CanSolveForVectorWithTallRandomMatrix()
        {
            var matrixA = Matrix<Complex32>.Build.Random(20, 10, 1);
            var matrixACopy = matrixA.Clone();
            var factorQR = matrixA.GramSchmidt();

            var vectorB = Vector<Complex32>.Build.Random(20, 1);
            var vectorX = factorQR.Solve(vectorB);

            // The solution x dimension is equal to the column dimension of A
            Assert.AreEqual(matrixA.ColumnCount, vectorX.Count);

            var test = (matrixA.ConjugateTranspose() * matrixA).Inverse() * matrixA.ConjugateTranspose() * vectorB;

            AssertHelpers.AlmostEqual(test, vectorX, 5);

            // Make sure A didn't change.
            AssertHelpers.AreEqual(matrixACopy, matrixA);
        }
    }
}
