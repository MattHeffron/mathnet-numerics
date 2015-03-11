// <copyright file="FailureStopCriterionTest.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.OneBased.Integer;
using MathNet.Numerics.LinearAlgebra.OneBased.Solvers;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Integer.Solvers.StopCriterion
{
    /// <summary>
    /// Failure stop criterion tests.
    /// </summary>
    [TestFixture, Category("LASolver")]
    public sealed class FailureStopCriterionTest
    {
        /// <summary>
        /// Can create.
        /// </summary>
        [Test]
        public void Create()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "Should have a criterion now");
        }

        /// <summary>
        /// Determine status with illegal iteration number throws <c>ArgumentOutOfRangeException</c>.
        /// </summary>
        [Test]
        public void DetermineStatusWithIllegalIterationNumberThrowsArgumentOutOfRangeException()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "There should be a criterion");

            Assert.That(() => criterion.DetermineStatus(-1, Vector<int>.Build.Dense(3, 4), Vector<int>.Build.Dense(3, 5), Vector<int>.Build.Dense(3, 6)), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// Determine status with non-matching vectors throws <c>ArgumentException</c>.
        /// </summary>
        [Test]
        public void DetermineStatusWithNonMatchingVectorsThrowsArgumentException()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "There should be a criterion");

            Assert.That(() => criterion.DetermineStatus(1, Vector<int>.Build.Dense(3, 4), Vector<int>.Build.Dense(3, 6), Vector<int>.Build.Dense(4, 4)), Throws.ArgumentException);
        }

		/////// <summary>
		/////// Can determine status with residual NaN.
		/////// </summary>
		////[Test]
		////public void DetermineStatusWithResidualNaN()
		////{
		////	var criterion = new FailureStopCriterion<int>();
		////	Assert.IsNotNull(criterion, "There should be a criterion");

		////	var solution = new DenseVector(new[] { 1, 1, 2 });
		////	var source = new DenseVector(new[] { 1001, 0, 2003 });
		////	var residual = new DenseVector(new[] { 1000, float.NaN, 2001 });

		////	var status = criterion.DetermineStatus(5, solution, source, residual);
		////	Assert.AreEqual(IterationStatus.Failure, status, "Should be failed");
		////}

		/////// <summary>
		/////// Can determine status with solution NaN.
		/////// </summary>
		////[Test]
		////public void DetermineStatusWithSolutionNaN()
		////{
		////	var criterion = new FailureStopCriterion<int>();
		////	Assert.IsNotNull(criterion, "There should be a criterion");

		////	var solution = new DenseVector(new[] { 1, 1, float.NaN });
		////	var source = new DenseVector(new[] { 1001, 0, 2003 });
		////	var residual = new DenseVector(new[] { 1000, 1000, 2001 });

		////	var status = criterion.DetermineStatus(5, solution, source, residual);
		////	Assert.AreEqual(IterationStatus.Failure, status, "Should be failed");
		////}

        /// <summary>
        /// Can determine status.
        /// </summary>
        [Test]
        public void DetermineStatus()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "There should be a criterion");

            var solution = new DenseVector(new[] { 3, 2, 1 });
            var source = new DenseVector(new[] { 1001, 0, 2003 });
            var residual = new DenseVector(new[] { 1, 2, 3 });

            var status = criterion.DetermineStatus(5, solution, source, residual);
            Assert.AreEqual(IterationStatus.Continue, status, "Should be running");
        }

        /// <summary>
        /// Can reset calculation state.
        /// </summary>
        [Test]
        public void ResetCalculationState()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "There should be a criterion");

            var solution = new DenseVector(new[] { 1, 1, 2 });
            var source = new DenseVector(new[] { 1001, 0, 2003 });
            var residual = new DenseVector(new[] { 1000, 1000, 2001 });

            var status = criterion.DetermineStatus(5, solution, source, residual);
            Assert.AreEqual(IterationStatus.Continue, status, "Should be running");

            criterion.Reset();
            Assert.AreEqual(IterationStatus.Continue, criterion.Status, "Should not have started");
        }

        /// <summary>
        /// Can clone stop criterion.
        /// </summary>
        [Test]
        public void Clone()
        {
            var criterion = new FailureStopCriterion<int>();
            Assert.IsNotNull(criterion, "There should be a criterion");

            var clone = criterion.Clone();
            Assert.IsInstanceOf(typeof (FailureStopCriterion<int>), clone, "Wrong criterion type");
        }
    }
}