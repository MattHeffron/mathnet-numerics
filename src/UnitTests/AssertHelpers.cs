// <copyright file="AssertHelpers.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
// Copyright (c) 2009-2010 Math.NET
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

using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests
{
#if NOSYSNUMERICS
    using Complex = Numerics.Complex;
#else
    using Complex = System.Numerics.Complex;
#endif

    /// <summary>
    /// A class which includes some assertion helper methods particularly for numerical code.
    /// </summary>
    internal static class AssertHelpers
    {
        public static void AlmostEqual(Complex expected, Complex actual)
        {
            if (expected.IsNaN() && actual.IsNaN() || expected.IsInfinity() && expected.IsInfinity())
            {
                return;
            }

            if (!expected.Real.AlmostEqual(actual.Real))
            {
                Assert.Fail("Real components are not equal. Expected:{0}; Actual:{1}", expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqual(actual.Imaginary))
            {
                Assert.Fail("Imaginary components are not equal. Expected:{0}; Actual:{1}", expected.Imaginary, actual.Imaginary);
            }
        }

        public static void AlmostEqual(Complex32 expected, Complex32 actual)
        {
            if (expected.IsNaN() && actual.IsNaN() || expected.IsInfinity() && expected.IsInfinity())
            {
                return;
            }

            if (!expected.Real.AlmostEqual(actual.Real))
            {
                Assert.Fail("Real components are not equal. Expected:{0}; Actual:{1}", expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqual(actual.Imaginary))
            {
                Assert.Fail("Imaginary components are not equal. Expected:{0}; Actual:{1}", expected.Imaginary, actual.Imaginary);
            }
        }


        public static void AlmostEqual(double expected, double actual, int decimalPlaces)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual))
            {
                return;
            }

            if (!expected.AlmostEqual(actual, decimalPlaces))
            {
                Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected, actual);
            }
        }

        public static void AlmostEqual(float expected, float actual, int decimalPlaces)
        {
            if (float.IsNaN(expected) && float.IsNaN(actual))
            {
                return;
            }

            if (!expected.AlmostEqual(actual, decimalPlaces))
            {
                Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected, actual);
            }
        }

        public static void AlmostEqual(Complex expected, Complex actual, int decimalPlaces)
        {
            if (!expected.Real.AlmostEqual(actual.Real, decimalPlaces))
            {
                Assert.Fail("Real components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqual(actual.Imaginary, decimalPlaces))
            {
                Assert.Fail("Imaginary components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Imaginary, actual.Imaginary);
            }
        }

        public static void AlmostEqual(Complex32 expected, Complex32 actual, int decimalPlaces)
        {
            if (!expected.Real.AlmostEqual(actual.Real, decimalPlaces))
            {
                Assert.Fail("Real components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqual(actual.Imaginary, decimalPlaces))
            {
                Assert.Fail("Imaginary components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Imaginary, actual.Imaginary);
            }
        }

        public static void AlmostEqualRelative(double expected, double actual, int decimalPlaces)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual))
            {
                return;
            }

            if (!expected.AlmostEqualRelative(actual, decimalPlaces))
            {
                Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected, actual);
            }
        }

        public static void AlmostEqualRelative(float expected, float actual, int decimalPlaces)
        {
            if (float.IsNaN(expected) && float.IsNaN(actual))
            {
                return;
            }

            if (!expected.AlmostEqualRelative(actual, decimalPlaces))
            {
                Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected, actual);
            }
        }

        public static void AlmostEqualRelative(Complex expected, Complex actual, int decimalPlaces)
        {
            if (!expected.Real.AlmostEqualRelative(actual.Real, decimalPlaces))
            {
                Assert.Fail("Real components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqualRelative(actual.Imaginary, decimalPlaces))
            {
                Assert.Fail("Imaginary components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Imaginary, actual.Imaginary);
            }
        }

        public static void AlmostEqualRelative(Complex32 expected, Complex32 actual, int decimalPlaces)
        {
            if (!expected.Real.AlmostEqualRelative(actual.Real, decimalPlaces))
            {
                Assert.Fail("Real components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Real, actual.Real);
            }

            if (!expected.Imaginary.AlmostEqualRelative(actual.Imaginary, decimalPlaces))
            {
                Assert.Fail("Imaginary components are not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.Imaginary, actual.Imaginary);
            }
        }


        public static void AlmostEqual(IList<double> expected, IList<double> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqual(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqual(IList<float> expected, IList<float> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqual(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqual(IList<Complex> expected, IList<Complex> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqual(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqual(IList<Complex32> expected, IList<Complex32> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqual(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        /// <remarks>
        /// There's no "Almost" about integer comparisons
        /// </remarks>
        public static void AreEqual(IList<int> expected, IList<int> actual)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                Assert.AreEqual(expected[i], actual[i]);
            }
        }

        public static void AlmostEqualRelative(IList<double> expected, IList<double> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqualRelative(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqualRelative(IList<float> expected, IList<float> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqualRelative(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqualRelative(IList<Complex> expected, IList<Complex> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqualRelative(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqualRelative(IList<Complex32> expected, IList<Complex32> actual, int decimalPlaces)
        {
            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual[i].AlmostEqualRelative(expected[i], decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected[i], actual[i]);
                }
            }
        }

        public static void AlmostEqual(Matrix<double> expected, Matrix<double> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(Matrix<float> expected, Matrix<float> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(Matrix<Complex> expected, Matrix<Complex> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(Matrix<Complex32> expected, Matrix<Complex32> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        /// <remarks>
        /// There's no "Almost" about integer comparisons
        /// </remarks>
        public static void AreEqual(Matrix<int> expected, Matrix<int> actual)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    Assert.AreEqual(expected.At(i, j), actual.At(i, j));
                }
            }
        }

        /// <remarks>
        /// There's no "Almost" about integer comparisons
        /// </remarks>
        public static void AreEqual(Vector<int> expected, Vector<int> actual)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                Assert.AreEqual(expected.At(i), actual.At(i));
            }
        }

        public static void AlmostEqual(Vector<double> expected, Vector<double> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(Vector<float> expected, Vector<float> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(Vector<Complex> expected, Vector<Complex> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(Vector<Complex32> expected, Vector<Complex32> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(Matrix<double> expected, Matrix<double> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(Matrix<float> expected, Matrix<float> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(Matrix<Complex> expected, Matrix<Complex> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(Matrix<Complex32> expected, Matrix<Complex32> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(Vector<double> expected, Vector<double> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(Vector<float> expected, Vector<float> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(Vector<Complex> expected, Vector<Complex> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(Vector<Complex32> expected, Vector<Complex32> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        #region OneBased
        #region Full Matrix and Vector checks
        // These could easily be copied and modified for use outside of OneBased...

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void AreEqual<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> actual) where T : struct, IEquatable<T>, IFormattable
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 1; i <= expected.RowCount; i++)
            {
                for (var j = 1; j <= expected.ColumnCount; j++)
                {
                    Assert.AreEqual(expected.At(i, j), actual.At(i, j));
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void AreEqual<T>(MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> actual) where T : struct, IEquatable<T>, IFormattable
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 1; i <= expected.Count; i++)
            {
                Assert.AreEqual(expected.At(i), actual.At(i));
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void IsUpperTriangular<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixU, T zero = default(T)) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything below the diagonal is zero
            for (var i = 1; i <= matrixU.RowCount; i++)
            {
                for (var j = 1; j < i; j++)
                {
                    Assert.AreEqual(zero, matrixU[i, j]);
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void IsLowerTriangular<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixL, T zero = default(T)) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything above the diagonal is zero
            for (var i = 1; i <= matrixL.RowCount; i++)
            {
                for (var j = i + 1; j <= matrixL.ColumnCount; j++)
                {
                    Assert.AreEqual(zero, matrixL[i, j]);
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void IsDiagonal<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixD, T zero = default(T)) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything off the diagonal is zero
            for (var i = 1; i <= matrixD.RowCount; i++)
            {
                for (var j = 1; j <= matrixD.ColumnCount; j++)
                {
                    if (i != j)
                    {
                        Assert.AreEqual(zero, matrixD[i, j]);
                    }
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void DiagonalHasValue<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixD, T value) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything on the diagonal is value
            int size = Math.Min(matrixD.RowCount, matrixD.ColumnCount);
            for (var i = 1; i <= size; i++)
            {
                Assert.AreEqual(value, matrixD[i, i]);
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void DiagonalValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixD, Action<T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Execute the valueAssertion on the diagonal
            int size = Math.Min(matrixD.RowCount, matrixD.ColumnCount);
            for (var i = 1; i <= size; i++)
            {
                valueAssertion(matrixD[i, i]);
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void DiagonalValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrixD, Action<int, T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Execute the valueAssertion on the diagonal
            int size = Math.Min(matrixD.RowCount, matrixD.ColumnCount);
            for (var i = 1; i <= size; i++)
            {
                valueAssertion(i, matrixD[i, i]);
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void ValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrix, Action<T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Execute the valueAssertion on all the elements of the matrix
            for (var i = 1; i <= matrix.RowCount; i++)
            {
                for (var j = 1; j <= matrix.ColumnCount; j++)
                {
                    valueAssertion(matrix[i, j]);
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void ValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<T> matrix, Action<int, int, T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Execute the valueAssertion on all the elements of the matrix
            for (var r = 1; r <= matrix.RowCount; r++)
            {
                for (var c = 1; c <= matrix.ColumnCount; c++)
                {
                    valueAssertion(r, c, matrix[r, c]);
                }
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void VectorHasValue<T>(MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> vector, T value) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything in the vector is value
            for (var i = 1; i <= vector.Count; i++)
            {
                Assert.AreEqual(value, vector[i]);
            }
        }

        /////// <remarks>
        /////// Simplify the unit tests to avoid needing to write the loop structures every time...
        /////// </remarks>
        ////public static void VectorValuesCondition<T>(MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> vector, Func<T, bool> valueTest) where T : struct, IEquatable<T>, IFormattable
        ////{
        ////    // Everything in the vector meets condition
        ////    for (var i = 1; i <= vector.Count; i++)
        ////    {
        ////        Assert.IsTrue(valueTest(vector[i]));
        ////    }
        ////}

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void ValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> vector, Action<T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything in the vector meets condition
            for (var i = 1; i <= vector.Count; i++)
            {
                valueAssertion(vector[i]);
            }
        }

        /// <remarks>
        /// Simplify the unit tests to avoid needing to write the loop structures every time...
        /// </remarks>
        public static void ValuesAssertion<T>(MathNet.Numerics.LinearAlgebra.OneBased.Vector<T> vector, Action<int, T> valueAssertion) where T : struct, IEquatable<T>, IFormattable
        {
            // Everything in the vector meets condition
            for (var i = 1; i <= vector.Count; i++)
            {
                valueAssertion(i, vector[i]);
            }
        }
        #endregion

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<double> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<double> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<float> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<float> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex32> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex32> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqual(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Vector<double> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<double> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Vector<float> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<float> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqual(MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex32> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex32> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqual(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<double> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<double> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<float> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<float> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex32> expected, MathNet.Numerics.LinearAlgebra.OneBased.Matrix<Complex32> actual, int decimalPlaces)
        {
            if (expected.ColumnCount != actual.ColumnCount || expected.RowCount != actual.RowCount)
            {
                Assert.Fail("Matrix dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.RowCount; i++)
            {
                for (var j = 0; j < expected.ColumnCount; j++)
                {
                    if (!actual.At(i, j).AlmostEqualRelative(expected.At(i, j), decimalPlaces))
                    {
                        Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i, j), actual.At(i, j));
                    }
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Vector<double> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<double> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Vector<float> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<float> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }

        public static void AlmostEqualRelative(MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex32> expected, MathNet.Numerics.LinearAlgebra.OneBased.Vector<Complex32> actual, int decimalPlaces)
        {
            if (expected.Count != actual.Count)
            {
                Assert.Fail("Vector dimensions mismatch. Expected: {0}; Actual: {1}", expected.ToTypeString(), actual.ToTypeString());
            }

            for (var i = 0; i < expected.Count; i++)
            {
                if (!actual.At(i).AlmostEqualRelative(expected.At(i), decimalPlaces))
                {
                    Assert.Fail("Not equal within {0} relative places. Expected:{1}; Actual:{2}", decimalPlaces, expected.At(i), actual.At(i));
                }
            }
        }
        #endregion
    }
}
