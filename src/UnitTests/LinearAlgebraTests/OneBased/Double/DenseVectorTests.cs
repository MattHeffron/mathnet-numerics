// <copyright file="DenseVectorTests.cs" company="Math.NET">
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

using MathNet.Numerics.LinearAlgebra.OneBased;
using MathNet.Numerics.LinearAlgebra.OneBased.Double;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Double
{
    /// <summary>
    /// Dense vector tests.
    /// </summary>
    public class DenseVectorTests : VectorTests
    {
        /// <summary>
        /// Creates a new instance of the Vector class.
        /// </summary>
        /// <param name="size">The size of the <strong>Vector</strong> to construct.</param>
        /// <returns>The new <c>Vector</c>.</returns>
        protected override Vector<double> CreateVector(int size)
        {
            return Vector<double>.Build.Dense(size);
        }

        /// <summary>
        /// Creates a new instance of the Vector class.
        /// </summary>
        /// <param name="data">The array to create this vector from.</param>
        /// <returns>The new <c>Vector</c>.</returns>
        protected override Vector<double> CreateVector(IList<double> data)
        {
            var vector = Vector<double>.Build.Dense(data.Count);
            for (var index = 0; index < data.Count; index++)
            {
                vector[index + 1] = data[index];
            }

            return vector;
        }

        /// <summary>
        /// Can create a dense vector form array.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromArray()
        {
            var data = new double[Data.Length];
            Array.Copy(Data, data, Data.Length);
            var vector = Vector<double>.Build.Dense(data);

            for (var i = 0; i < data.Length; i++)
            {
                Assert.AreEqual(data[i], vector[i + 1]);
            }

            vector[1] = 100.0;
            Assert.AreEqual(100.0, data[0]);
        }

        /// <summary>
        /// Can create a dense vector from another dense vector.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromAnotherDenseVector()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = Vector<double>.Build.DenseOfVector(vector);

            Assert.AreNotSame(vector, other);
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a dense vector from another vector.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromAnotherVector()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = Vector<double>.Build.DenseOfVector(vector);

            Assert.AreNotSame(vector, other);
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a dense vector from user defined vector.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromUserDefinedVector()
        {
            var vector = new UserDefinedVector(Data);
            var other = Vector<double>.Build.DenseOfVector(vector);

            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a dense vector with constant values.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorWithConstantValues()
        {
            var vector = Vector<double>.Build.Dense(5, 5);
            AssertHelpers.AllVectorElementsHaveValue(vector, 5);
        }

        /// <summary>
        /// Can create a dense matrix.
        /// </summary>
        [Test]
        public void CanCreateDenseMatrix()
        {
            var vector = Vector<double>.Build.Dense(3);
            var matrix = Matrix<double>.Build.SameAs(vector, 2, 3);
            Assert.IsInstanceOf<DenseMatrix>(matrix);
            Assert.AreEqual(2, matrix.RowCount);
            Assert.AreEqual(3, matrix.ColumnCount);
        }

        /// <summary>
        /// Can convert a dense vector to an array.
        /// </summary>
        [Test]
        public void CanConvertDenseVectorToArray()
        {
            var vector = new DenseVector(Data);
            var array = (double[])vector;
            Assert.IsInstanceOf(typeof(double[]), array);
            CollectionAssert.AreEqual(vector, array);
        }

        /// <summary>
        /// Can convert an array to a dense vector.
        /// </summary>
        [Test]
        public void CanConvertArrayToDenseVector()
        {
            var array = new[] { 0.0, 1.0, 2.0, 3.0, 4.0 };
            var vector = (DenseVector)array;
            Assert.IsInstanceOf(typeof(DenseVector), vector);
            CollectionAssert.AreEqual(array, array);
        }

        /// <summary>
        /// Can call unary plus operator on a vector.
        /// </summary>
        [Test]
        public void CanCallUnaryPlusOperatorOnDenseVector()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = +vector;
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can add two dense vectors using "+" operator.
        /// </summary>
        [Test]
        public void CanAddTwoDenseVectorsUsingOperator()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = Vector<double>.Build.Dense(Data);
            var result = vector + other;
            CollectionAssert.AreEqual(Data, vector, "Making sure the original vector wasn't modified.");
            CollectionAssert.AreEqual(Data, other, "Making sure the original vector wasn't modified.");
            AssertHelpers.IndexedAssertion(result, i => Assert.AreEqual(Data[i - 1]*2.0, result[i]));
        }

        /// <summary>
        /// Can call unary negate operator on a dense vector.
        /// </summary>
        [Test]
        public void CanCallUnaryNegationOperatorOnDenseVector()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = -vector;
            AssertHelpers.IndexedAssertion(other, i => Assert.AreEqual(-Data[i - 1], other[i]));
        }

        /// <summary>
        /// Can subtract two dense vectors using "-" operator.
        /// </summary>
        [Test]
        public void CanSubtractTwoDenseVectorsUsingOperator()
        {
            var vector = Vector<double>.Build.Dense(Data);
            var other = Vector<double>.Build.Dense(Data);
            var result = vector - other;
            CollectionAssert.AreEqual(Data, vector, "Making sure the original vector wasn't modified.");
            CollectionAssert.AreEqual(Data, other, "Making sure the original vector wasn't modified.");
            AssertHelpers.AllVectorElementsHaveValue(result, 0.0);
        }

        /// <summary>
        /// Can multiply a dense vector by a scalar using "*" operator.
        /// </summary>
        [Test]
        public void CanMultiplyDenseVectorByScalarUsingOperators()
        {
            var vector = Vector<double>.Build.Dense(Data);
            vector = vector * 2.0;
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] * 2.0, vector[i]));

            vector = vector * 1.0;
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] * 2.0, vector[i]));

            vector = Vector<double>.Build.Dense(Data);
            vector = 2.0 * vector;

            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] * 2.0, vector[i]));

            vector = 1.0 * vector;
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] * 2.0, vector[i]));
        }

        /// <summary>
        /// Can divide a dense vector by a scalar using "/" operator.
        /// </summary>
        [Test]
        public void CanDivideDenseVectorByScalarUsingOperators()
        {
            var vector = Vector<double>.Build.Dense(Data);
            vector = vector / 2.0;
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] / 2.0, vector[i]));
            for (var i = 0; i < Data.Length; i++)
            {
                Assert.AreEqual(Data[i] / 2.0, vector[i]);
            }

            vector = vector / 1.0;
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(Data[i - 1] / 2.0, vector[i]));

        }

        /// <summary>
        /// Can calculate an outer product for a dense vector.
        /// </summary>
        [Test]
        public void CanCalculateOuterProductForDenseVector()
        {
            var vector1 = CreateVector(Data);
            var vector2 = CreateVector(Data);
            var m = Vector<double>.OuterProduct(vector1, vector2);
            AssertHelpers.IndexedAssertion(m, (i, j) => Assert.AreEqual(vector1[i] * vector2[j], m[i, j]));
        }

        [Test]
        public void VectorToVectorString()
        {
            var v = Vector<double>.Build.Dense(20);
            for (int i = 1; i < 25; i++)
            {
                GC.KeepAlive(v.ToVectorString(i, 80));
            }
        }
    }
}
