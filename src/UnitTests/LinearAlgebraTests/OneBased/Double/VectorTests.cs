// <copyright file="VectorTests.cs" company="Math.NET">
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
using System.Collections;
using System.Collections.Generic;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra.OneBased;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Double
{
    /// <summary>
    /// Abstract class with the common set of vector tests.
    /// </summary>
    public abstract partial class VectorTests
    {
        /// <summary>
        /// Test vector values.
        /// </summary>
        protected readonly double[] Data = { 1, 2, 3, 4, 5 };

        /// <summary>
        /// Can clone a vector.
        /// </summary>
        [Test]
        public void CanCloneVector()
        {
            var vector = CreateVector(Data);
            var clone = vector.Clone();

            Assert.AreNotSame(vector, clone);
            Assert.AreEqual(vector.Count, clone.Count);
            AssertHelpers.AreEqual(vector, clone);
        }

#if !PORTABLE
        /// <summary>
        /// Can clone a vector using <c>IClonable</c> interface method.
        /// </summary>
        [Test]
        public void CanCloneVectorUsingICloneable()
        {
            var vector = CreateVector(Data);
            var clone = (Vector<double>)((ICloneable)vector).Clone();

            Assert.AreNotSame(vector, clone);
            Assert.AreEqual(vector.Count, clone.Count);
            AssertHelpers.AreEqual(vector, clone);
        }
#endif

        /// <summary>
        /// Can copy part of a vector to another vector.
        /// </summary>
        [Test]
        public void CanCopyPartialVectorToAnother()
        {
            var vector = CreateVector(Data);
            var other = CreateVector(Data.Length);

            vector.CopySubVectorTo(other, 3, 3, 2);

            Assert.AreEqual(0.0, other[1]);
            Assert.AreEqual(0.0, other[2]);
            Assert.AreEqual(3.0, other[3]);
            Assert.AreEqual(4.0, other[4]);
            Assert.AreEqual(0.0, other[5]);
        }

        /// <summary>
        /// Can copy part of a vector to the same vector.
        /// </summary>
        [Test]
        public void CanCopyPartialVectorToSelf()
        {
            var vector = CreateVector(Data);
            vector.CopySubVectorTo(vector, 1, 3, 2);

            Assert.AreEqual(1.0, vector[1]);
            Assert.AreEqual(2.0, vector[2]);
            Assert.AreEqual(1.0, vector[3]);
            Assert.AreEqual(2.0, vector[4]);
            Assert.AreEqual(5.0, vector[5]);
        }

        /// <summary>
        /// Can copy a vector to another vector.
        /// </summary>
        [Test]
        public void CanCopyVectorToAnother()
        {
            var vector = CreateVector(Data);
            var other = CreateVector(Data.Length);

            vector.CopyTo(other);
            AssertHelpers.AreEqual(vector, other);
        }

        /// <summary>
        /// Can create a matrix using instance of a vector.
        /// </summary>
        [Test]
        public void CanCreateMatrix()
        {
            var vector = CreateVector(Data);
            var matrix = Matrix<double>.Build.SameAs(vector, 10, 10);
            Assert.AreEqual(matrix.RowCount, 10);
            Assert.AreEqual(matrix.ColumnCount, 10);
        }

        /// <summary>
        /// Can create a vector using the instance of vector.
        /// </summary>
        [Test]
        public void CanCreateVector()
        {
            var expected = CreateVector(5);
            var actual = Vector<double>.Build.SameAs(expected, 5);
            Assert.AreEqual(expected.Storage.IsDense, actual.Storage.IsDense, "vectors are same kind.");
        }

        /// <summary>
        /// Can enumerate over a vector.
        /// </summary>
        [Test]
        public void CanEnumerateOverVector()
        {
            var vector = CreateVector(Data);
            for (var i = 0; i < Data.Length; i++)
            {
                Assert.AreEqual(Data[i], vector[i + 1]);
            }
        }

        /// <summary>
        /// Can enumerate over a vector using <c>IEnumerable</c> interface.
        /// </summary>
        [Test]
        public void CanEnumerateOverVectorUsingIEnumerable()
        {
            var enumerable = (IEnumerable)CreateVector(Data);
            var index = 0;
            foreach (var element in enumerable)
            {
                Assert.AreEqual(Data[index++], (double)element);
            }
        }

        /// <summary>
        /// Can equate vectors.
        /// </summary>
        [Test]
        public void CanEquateVectors()
        {
            var vector1 = CreateVector(Data);
            var vector2 = CreateVector(Data);
            var vector3 = CreateVector(4);
            Assert.IsTrue(vector1.Equals(vector1));
            Assert.IsTrue(vector1.Equals(vector2));
            Assert.IsFalse(vector1.Equals(vector3));
            Assert.IsFalse(vector1.Equals(null));
            Assert.IsFalse(vector1.Equals(2));
        }

        /// <summary>
        /// <c>CreateVector</c> throws <c>ArgumentOutOfRangeException</c> if size is negative.
        /// </summary>
        [Test]
        public void SizeIsNegativeThrowsArgumentOutOfRangeException()
        {
            Assert.That(() => CreateVector(-1), Throws.TypeOf<ArgumentOutOfRangeException>());
            // 0 length is allowed on OneBased because Matlab allows it
            //Assert.That(() => CreateVector(0), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// Can equate vectors using Object.Equals.
        /// </summary>
        [Test]
        public void CanEquateVectorsUsingObjectEquals()
        {
            var vector1 = CreateVector(Data);
            var vector2 = CreateVector(Data);
            Assert.IsTrue(vector1.Equals((object)vector2));
        }

        /// <summary>
        /// Can get hash code of a vector.
        /// </summary>
        [Test]
        public void CanGetHashCode()
        {
            var vector = CreateVector(new double[] { 1, 2, 3, 4 });
            Assert.AreEqual(vector.GetHashCode(), vector.GetHashCode());
            Assert.AreEqual(vector.GetHashCode(), CreateVector(new double[] { 1, 2, 3, 4 }).GetHashCode());
            Assert.AreNotEqual(vector.GetHashCode(), CreateVector(new double[] { 1 }).GetHashCode());
        }

        /// <summary>
        /// Can enumerate over a vector using indexed enumerator.
        /// </summary>
        [Test]
        public void CanEnumerateOverVectorUsingIndexedEnumerator()
        {
            var vector = CreateVector(Data);
            foreach (var pair in vector.EnumerateIndexed())
            {
                Assert.AreEqual(Data[pair.Item1 - 1], pair.Item2);
            }
        }

        /// <summary>
        /// Can enumerate over a vector using non-zero enumerator.
        /// </summary>
        [Test]
        public void CanEnumerateOverVectorUsingNonZeroEnumerator()
        {
            var vector = CreateVector(Data);
            foreach (var pair in vector.EnumerateIndexed(Zeros.AllowSkip))
            {
                Assert.AreEqual(Data[pair.Item1 - 1], pair.Item2);
                Assert.AreNotEqual(0d, pair.Item2);
            }
        }

        /// <summary>
        /// Can convert a vector to array.
        /// </summary>
        [Test]
        public void CanConvertVectorToArray()
        {
            var vector = CreateVector(Data);
            var array = vector.ToArray();
            CollectionAssert.AreEqual(array, vector);
        }

        /// <summary>
        /// Can convert a vector to column matrix.
        /// </summary>
        [Test]
        public void CanConvertVectorToColumnMatrix()
        {
            var vector = CreateVector(Data);
            var matrix = vector.ToColumnMatrix();

            Assert.AreEqual(vector.Count, matrix.RowCount);
            Assert.AreEqual(1, matrix.ColumnCount);
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(vector[i], matrix[i, 1]));
        }

        /// <summary>
        /// Can convert a vector to row matrix.
        /// </summary>
        [Test]
        public void CanConvertVectorToRowMatrix()
        {
            var vector = CreateVector(Data);
            var matrix = vector.ToRowMatrix();

            Assert.AreEqual(vector.Count, matrix.ColumnCount);
            Assert.AreEqual(1, matrix.RowCount);
            AssertHelpers.IndexedAssertion(vector, i => Assert.AreEqual(vector[i], matrix[1, i]));
        }

        /// <summary>
        /// Can set values in vector.
        /// </summary>
        [Test]
        public void CanSetValues()
        {
            var vector = CreateVector(Data);
            vector.SetValues(Data);
            CollectionAssert.AreEqual(vector, Data);
        }

        /// <summary>
        /// Can get subvector from a vector.
        /// </summary>
        /// <param name="index">The first element to begin copying from.</param>
        /// <param name="length">The number of elements to copy.</param>
        [TestCase(1, 5)]
        [TestCase(3, 2)]
        [TestCase(2, 4)]
        public void CanGetSubVector(int index, int length)
        {
            var vector = CreateVector(Data);
            var sub = vector.SubVector(index, length);
            Assert.AreEqual(length, sub.Count);
            AssertHelpers.IndexedAssertion(sub, i => Assert.AreEqual(vector[i + index - 1], sub[i]));
        }

        /// <summary>
        /// Getting subvector using wrong parameters throw an exception.
        /// </summary>
        /// <param name="index">The first element to begin copying from.</param>
        /// <param name="length">The number of elements to copy.</param>
        [TestCase(6, 10)]
        [TestCase(1, -10)]
        public void CanGetSubVectorWithWrongValuesShouldThrowException(int index, int length)
        {
            var vector = CreateVector(Data);
            Assert.That(() => vector.SubVector(index, length), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// Can find absolute minimum value index.
        /// </summary>
        [Test]
        public void CanFindAbsoluteMinimumIndex()
        {
            var source = CreateVector(Data);
            const int Expected = 1;
            var actual = source.AbsoluteMinimumIndex();
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find absolute minimum value of a vector.
        /// </summary>
        [Test]
        public void CanFindAbsoluteMinimum()
        {
            var source = CreateVector(Data);
            const double Expected = 1;
            var actual = source.AbsoluteMinimum();
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find absolute maximum value index.
        /// </summary>
        [Test]
        public void CanFindAbsoluteMaximumIndex()
        {
            var source = CreateVector(Data);
            const int Expected = 5;
            var actual = source.AbsoluteMaximumIndex();
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find absolute maximum value of a vector.
        /// </summary>
        [Test]
        public void CanFindAbsoluteMaximum()
        {
            var source = CreateVector(Data);
            const double Expected = 5;
            var actual = source.AbsoluteMaximum();
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find maximum value index.
        /// </summary>
        [Test]
        public void CanFindMaximumIndex()
        {
            var vector = CreateVector(Data);

            const int Expected = 5;
            var actual = vector.MaximumIndex();

            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find maximum value of a vector.
        /// </summary>
        [Test]
        public void CanFindMaximum()
        {
            var vector = CreateVector(Data);

            const double Expected = 5;
            var actual = vector.Maximum();

            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find minimum value index.
        /// </summary>
        [Test]
        public void CanFindMinimumIndex()
        {
            var vector = CreateVector(Data);

            const int Expected = 1;
            var actual = vector.MinimumIndex();

            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can find minimum value of a vector.
        /// </summary>
        [Test]
        public void CanFindMinimum()
        {
            var vector = CreateVector(Data);

            const double Expected = 1;
            var actual = vector.Minimum();

            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can compute the sum of a vector elements.
        /// </summary>
        [Test]
        public void CanSum()
        {
            double[] testData = { -20, -10, 10, 20, 30 };
            var vector = CreateVector(testData);
            var actual = vector.Sum();
            const double Expected = 30;
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Can compute the sum of the absolute value a vector elements.
        /// </summary>
        [Test]
        public void CanSumMagnitudes()
        {
            double[] testData = { -20, -10, 10, 20, 30 };
            var vector = CreateVector(testData);
            var actual = vector.SumMagnitudes();
            const double Expected = 90;
            Assert.AreEqual(Expected, actual);
        }

        /// <summary>
        /// Set values with <c>null</c> parameter throw exception.
        /// </summary>
        [Test]
        public void SetValuesWithNullParameterThrowsArgumentException()
        {
            var vector = CreateVector(Data);
            Assert.That(() => vector.SetValues(null), Throws.TypeOf<ArgumentNullException>());
        }

        /// <summary>
        /// Set values with non-equal data length throw exception.
        /// </summary>
        [Test]
        public void SetValuesWithNonEqualDataLengthThrowsArgumentException()
        {
            var vector = CreateVector(Data.Length + 2);
            Assert.That(() => vector.SetValues(Data), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// Generate a vector with number of elements less than zero throw an exception.
        /// </summary>
        [Test]
        public void RandomWithNumberOfElementsLessThanZeroThrowsArgumentException()
        {
            Assert.That(() => Vector<double>.Build.Random(-2, new ContinuousUniform()), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// Can clear a vector.
        /// </summary>
        [Test]
        public void CanClearVector()
        {
            double[] testData = { -20, -10, 10, 20, 30 };
            var vector = CreateVector(testData);
            vector.Clear();
            AssertHelpers.AllVectorElementsHaveValue(vector, 0.0);
        }

        /// <summary>
        /// Creates a new instance of the Vector class.
        /// </summary>
        /// <param name="size">The size of the <strong>Vector</strong> to construct.</param>
        /// <returns>The new <c>Vector</c>.</returns>
        protected abstract Vector<double> CreateVector(int size);

        /// <summary>
        /// Creates a new instance of the Vector class.
        /// </summary>
        /// <param name="data">The array to create this vector from.</param>
        /// <returns>The new <c>Vector</c>.</returns>
        protected abstract Vector<double> CreateVector(IList<double> data);
    }
}
