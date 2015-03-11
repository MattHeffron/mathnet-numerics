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
using MathNet.Numerics.LinearAlgebra.OneBased.Complex;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Complex
{
#if NOSYSNUMERICS
    using Complex = Numerics.Complex;
#else
    using Complex = System.Numerics.Complex;
#endif

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
        protected override Vector<Complex> CreateVector(int size)
        {
            return new DenseVector(size);
        }

        /// <summary>
        /// Creates a new instance of the Vector class.
        /// </summary>
        /// <param name="data">The array to create this vector from.</param>
        /// <returns>The new <c>Vector</c>.</returns>
        protected override Vector<Complex> CreateVector(IList<Complex> data)
        {
            var vector = new DenseVector(data.Count);
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
            var data = new Complex[Data.Length];
            Array.Copy(Data, data, Data.Length);
            var vector = new DenseVector(data);

            for (var i = 0; i < data.Length; i++)
            {
                Assert.AreEqual(data[i], vector[i + 1]);
            }

            vector[1] = new Complex(10.0, 1);
            Assert.AreEqual(new Complex(10.0, 1), data[0]);
        }

        /// <summary>
        /// Can create a dense vector from another dense vector.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromAnotherDenseVector()
        {
            var vector = new DenseVector(Data);
            var other = DenseVector.OfVector(vector);

            Assert.AreNotSame(vector, other);
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a dense vector from another vector.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorFromAnotherVector()
        {
            var vector = (Vector<Complex>) new DenseVector(Data);
            var other = DenseVector.OfVector(vector);

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
            var other = DenseVector.OfVector(vector);

            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a dense vector with constant values.
        /// </summary>
        [Test]
        public void CanCreateDenseVectorWithConstantValues()
        {
            var vector = DenseVector.Create(5, 5);
            AssertHelpers.VectorHasValue(vector, new Complex(5.0, 0));
        }

        /// <summary>
        /// Can create a dense matrix.
        /// </summary>
        [Test]
        public void CanCreateDenseMatrix()
        {
            var vector = new DenseVector(3);
            var matrix = Matrix<Complex>.Build.SameAs(vector, 2, 3);
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
            var array = (Complex[]) vector;
            Assert.IsInstanceOf(typeof (Complex[]), array);
            CollectionAssert.AreEqual(vector, array);
        }

        /// <summary>
        /// Can convert an array to a dense vector.
        /// </summary>
        [Test]
        public void CanConvertArrayToDenseVector()
        {
            var array = new[] {new Complex(1, 1), new Complex(2, 1), new Complex(3, 1), new Complex(4, 1)};
            var vector = (DenseVector) array;
            Assert.IsInstanceOf(typeof (DenseVector), vector);
            CollectionAssert.AreEqual(array, array);
        }

        /// <summary>
        /// Can call unary plus operator on a vector.
        /// </summary>
        [Test]
        public void CanCallUnaryPlusOperatorOnDenseVector()
        {
            var vector = new DenseVector(Data);
            var other = +vector;
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can add two dense vectors using "+" operator.
        /// </summary>
        [Test]
        public void CanAddTwoDenseVectorsUsingOperator()
        {
            var vector = new DenseVector(Data);
            var other = new DenseVector(Data);
            var result = vector + other;
            CollectionAssert.AreEqual(Data, vector, "Making sure the original vector wasn't modified.");
            CollectionAssert.AreEqual(Data, other, "Making sure the original vector wasn't modified.");
            AssertHelpers.ValuesAssertion(result, (i, v) => Assert.AreEqual(Data[i - 1]*2.0, v));
        }

        /// <summary>
        /// Can call unary negate operator on a dense vector.
        /// </summary>
        [Test]
        public void CanCallUnaryNegationOperatorOnDenseVector()
        {
            var vector = new DenseVector(Data);
            var other = -vector;
            AssertHelpers.ValuesAssertion(other, (i, v) => Assert.AreEqual(-Data[i - 1], v));
        }

        /// <summary>
        /// Can subtract two dense vectors using "-" operator.
        /// </summary>
        [Test]
        public void CanSubtractTwoDenseVectorsUsingOperator()
        {
            var vector = new DenseVector(Data);
            var other = new DenseVector(Data);
            var result = vector - other;
            CollectionAssert.AreEqual(Data, vector, "Making sure the original vector wasn't modified.");
            CollectionAssert.AreEqual(Data, other, "Making sure the original vector wasn't modified.");
            AssertHelpers.VectorHasValue(result, Complex.Zero);
        }

        /// <summary>
        /// Can multiply a dense vector by a scalar using "*" operator.
        /// </summary>
        [Test]
        public void CanMultiplyDenseVectorByScalarUsingOperators()
        {
            var vector = new DenseVector(Data);
            vector = vector*new Complex(2.0, 1);
            AssertHelpers.ValuesAssertion(vector, (i, v) => Assert.AreEqual(Data[i - 1] * new Complex(2.0, 1), v));

            vector = vector*1.0;
            AssertHelpers.ValuesAssertion(vector, (i, v) => Assert.AreEqual(Data[i - 1] * new Complex(2.0, 1), v));

            vector = new DenseVector(Data);
            vector = new Complex(2.0, 1)*vector;
            AssertHelpers.ValuesAssertion(vector, (i, v) => Assert.AreEqual(Data[i - 1] * new Complex(2.0, 1), v));

            vector = 1.0*vector;
            AssertHelpers.ValuesAssertion(vector, (i, v) => Assert.AreEqual(Data[i - 1] * new Complex(2.0, 1), v));
        }

        /// <summary>
        /// Can divide a dense vector by a scalar using "/" operator.
        /// </summary>
        [Test]
        public void CanDivideDenseVectorByComplexUsingOperators()
        {
            var vector = new DenseVector(Data);
            vector = vector/new Complex(2.0, 1);
            AssertHelpers.ValuesAssertion(vector, (i, v) => AssertHelpers.AlmostEqualRelative(Data[i - 1] / new Complex(2.0, 1), v, 14));

            vector = vector/1.0;
            AssertHelpers.ValuesAssertion(vector, (i, v) => AssertHelpers.AlmostEqualRelative(Data[i - 1] / new Complex(2.0, 1), v, 14));
        }

        /// <summary>
        /// Can calculate an outer product for a dense vector.
        /// </summary>
        [Test]
        public void CanCalculateOuterProductForDenseVector()
        {
            var vector1 = CreateVector(Data);
            var vector2 = CreateVector(Data);
            var m = Vector<Complex>.OuterProduct(vector1, vector2);
            AssertHelpers.ValuesAssertion(m, (i, j, v) => Assert.AreEqual(vector1[i] * vector2[j], m[i, j]));
        }
    }
}