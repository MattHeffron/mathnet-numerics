// <copyright file="MlkBiCgStabTest.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.OneBased;
using MathNet.Numerics.LinearAlgebra.OneBased.Complex32;
using MathNet.Numerics.LinearAlgebra.OneBased.Complex32.Solvers;
using MathNet.Numerics.LinearAlgebra.OneBased.Solvers;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Complex32.Solvers.Iterative
{
    using Numerics;

    /// <summary>
    /// Tests for Multiple-Lanczos Bi-Conjugate Gradient stabilized iterative matrix solver.
    /// </summary>
    [TestFixture, Category("LASolver")]
    public class MlkBiCgStabTest
    {
        /// <summary>
        /// Convergence boundary.
        /// </summary>
        const float ConvergenceBoundary = 1e-5f;

        /// <summary>
        /// Maximum iterations.
        /// </summary>
        const int MaximumIterations = 1000;

        /// <summary>
        /// Solve wide matrix throws <c>ArgumentException</c>.
        /// </summary>
        [Test]
        public void SolveWideMatrixThrowsArgumentException()
        {
            var matrix = new SparseMatrix(2, 3);
            var input = new DenseVector(2);

            var solver = new MlkBiCgStab();
            Assert.That(() => matrix.SolveIterative(input, solver), Throws.ArgumentException);
        }

        /// <summary>
        /// Solve long matrix throws <c>ArgumentException</c>.
        /// </summary>
        [Test]
        public void SolveLongMatrixThrowsArgumentException()
        {
            var matrix = new SparseMatrix(3, 2);
            var input = new DenseVector(3);

            var solver = new MlkBiCgStab();
            Assert.That(() => matrix.SolveIterative(input, solver), Throws.ArgumentException);
        }

        /// <summary>
        /// Solve unit matrix and back multiply.
        /// </summary>
        [Test]
        public void SolveUnitMatrixAndBackMultiply()
        {
            // Create the identity matrix
            var matrix = SparseMatrix.CreateIdentity(100);

            // Create the y vector
            var y = Vector<Complex32>.Build.Dense(matrix.RowCount, 1);

            // Create an iteration monitor which will keep track of iterative convergence
            var monitor = new Iterator<Complex32>(
                new IterationCountStopCriterion<Complex32>(MaximumIterations),
                new ResidualStopCriterion<Complex32>(ConvergenceBoundary),
                new DivergenceStopCriterion<Complex32>(),
                new FailureStopCriterion<Complex32>());

            var solver = new MlkBiCgStab();

            // Solve equation Ax = y
            var x = matrix.SolveIterative(y, solver, monitor);

            // Now compare the results
            Assert.IsNotNull(x, "#02");
            Assert.AreEqual(y.Count, x.Count, "#03");

            // Back multiply the vector
            var z = matrix.Multiply(x);

            // Check that the solution converged
            Assert.IsTrue(monitor.Status == IterationStatus.Converged, "#04");

            // Now compare the vectors
            AssertHelpers.ValuesAssertion(y, (i, v) => Assert.GreaterOrEqual(ConvergenceBoundary, (y[i] - z[i]).Magnitude, "#05-" + i));
        }

        /// <summary>
        /// Solve scaled unit matrix and back multiply.
        /// </summary>
        [Test]
        public void SolveScaledUnitMatrixAndBackMultiply()
        {
            // Create the identity matrix
            var matrix = SparseMatrix.CreateIdentity(100);

            // Scale it with a funny number
            matrix.Multiply((float)Math.PI, matrix);

            // Create the y vector
            var y = Vector<Complex32>.Build.Dense(matrix.RowCount, 1);

            // Create an iteration monitor which will keep track of iterative convergence
            var monitor = new Iterator<Complex32>(
                new IterationCountStopCriterion<Complex32>(MaximumIterations),
                new ResidualStopCriterion<Complex32>(ConvergenceBoundary),
                new DivergenceStopCriterion<Complex32>(),
                new FailureStopCriterion<Complex32>());

            var solver = new MlkBiCgStab();

            // Solve equation Ax = y
            var x = matrix.SolveIterative(y, solver, monitor);

            // Now compare the results
            Assert.IsNotNull(x, "#02");
            Assert.AreEqual(y.Count, x.Count, "#03");

            // Back multiply the vector
            var z = matrix.Multiply(x);

            // Check that the solution converged
            Assert.IsTrue(monitor.Status == IterationStatus.Converged, "#04");

            // Now compare the vectors
            AssertHelpers.ValuesAssertion(y, (i, v) => Assert.GreaterOrEqual(ConvergenceBoundary, (y[i] - z[i]).Magnitude, "#05-" + i));
        }

        /// <summary>
        /// Solve poisson matrix and back multiply.
        /// </summary>
        [Test]
        public void SolvePoissonMatrixAndBackMultiply()
        {
            // Create the matrix
            var matrix = MatrixHelpers.MakePoissonTestMatrix<Complex32>();

            // Create the y vector
            var y = Vector<Complex32>.Build.Dense(matrix.RowCount, 1);

            // Create an iteration monitor which will keep track of iterative convergence
            var monitor = new Iterator<Complex32>(
                new IterationCountStopCriterion<Complex32>(MaximumIterations),
                new ResidualStopCriterion<Complex32>(ConvergenceBoundary),
                new DivergenceStopCriterion<Complex32>(),
                new FailureStopCriterion<Complex32>());

            var solver = new MlkBiCgStab();

            // Solve equation Ax = y
            var x = matrix.SolveIterative(y, solver, monitor);

            // Now compare the results
            Assert.IsNotNull(x, "#02");
            Assert.AreEqual(y.Count, x.Count, "#03");

            // Back multiply the vector
            var z = matrix.Multiply(x);

            // Check that the solution converged
            Assert.IsTrue(monitor.Status == IterationStatus.Converged, "#04");

            // Now compare the vectors
            AssertHelpers.ValuesAssertion(y, (i, v) => Assert.GreaterOrEqual(ConvergenceBoundary, (y[i] - z[i]).Magnitude, "#05-" + i));
        }

        /// <summary>
        /// Can solve for a random vector.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(4)]
        public void CanSolveForRandomVector(int order)
        {
            for (var iteration = 5; iteration > 3; iteration--)
            {
                var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
                var vectorb = Vector<Complex32>.Build.Random(order, 1);

                var monitor = new Iterator<Complex32>(
                    new IterationCountStopCriterion<Complex32>(1000),
                    new ResidualStopCriterion<Complex32>(Math.Pow(1.0/10.0, iteration)));

                var solver = new MlkBiCgStab();

                var resultx = matrixA.SolveIterative(vectorb, solver, monitor);

                if (monitor.Status != IterationStatus.Converged)
                {
                    // Solution was not found, try again downgrading convergence boundary
                    continue;
                }

                Assert.AreEqual(matrixA.ColumnCount, resultx.Count);
                var matrixBReconstruct = matrixA*resultx;

                // Check the reconstruction.
                AssertHelpers.AlmostEqual(vectorb, matrixBReconstruct, iteration - 3);

                return;
            }

            Assert.Fail("Solution was not found in 3 tries");
        }

        /// <summary>
        /// Can solve for random matrix.
        /// </summary>
        /// <param name="order">Matrix order.</param>
        [TestCase(4)]
        public void CanSolveForRandomMatrix(int order)
        {
            for (var iteration = 5; iteration > 3; iteration--)
            {
                var matrixA = Matrix<Complex32>.Build.Random(order, order, 1);
                var matrixB = Matrix<Complex32>.Build.Random(order, order, 1);

                var monitor = new Iterator<Complex32>(
                    new IterationCountStopCriterion<Complex32>(1000),
                    new ResidualStopCriterion<Complex32>(Math.Pow(1.0/10.0, iteration)));

                var solver = new MlkBiCgStab();
                var matrixX = matrixA.SolveIterative(matrixB, solver, monitor);

                if (monitor.Status != IterationStatus.Converged)
                {
                    // Solution was not found, try again downgrading convergence boundary
                    continue;
                }

                // The solution X row dimension is equal to the column dimension of A
                Assert.AreEqual(matrixA.ColumnCount, matrixX.RowCount);

                // The solution X has the same number of columns as B
                Assert.AreEqual(matrixB.ColumnCount, matrixX.ColumnCount);

                var matrixBReconstruct = matrixA*matrixX;

                // Check the reconstruction.
                AssertHelpers.AlmostEqual(matrixB, matrixBReconstruct, iteration - 3);

                return;
            }

            Assert.Fail("Solution was not found in 3 tries");
        }
    }
}
