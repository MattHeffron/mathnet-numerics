// <copyright file="ManagedLinearAlgebraProvider.Single.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Threading;

namespace MathNet.Numerics.Providers.LinearAlgebra
{

#if !NOSYSNUMERICS
    using Complex = System.Numerics.Complex;
#endif

    /// <summary>
    /// The managed linear algebra provider.
    /// </summary>
    public partial class ManagedLinearAlgebraProvider
    {
        /// <summary>
        /// Adds a scaled vector to another: <c>result = y + alpha*x</c>.
        /// </summary>
        /// <param name="y">The vector to update.</param>
        /// <param name="alpha">The value to scale <paramref name="x"/> by.</param>
        /// <param name="x">The vector to add to <paramref name="y"/>.</param>
        /// <param name="result">The result of the addition.</param>
        /// <remarks>This is similar to the AXPY BLAS routine.</remarks>
        public virtual void AddVectorToScaledVector(int[] y, int alpha, int[] x, int[] result)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (y.Length != x.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            if (alpha == 0)
            {
                y.Copy(result);
            }
            else if (alpha == 1)
            {
                CommonParallel.For(0, y.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        result[i] = y[i] + x[i];
                    }
                });
            }
            else
            {
                CommonParallel.For(0, y.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        result[i] = y[i] + (alpha * x[i]);
                    }
                });
            }
        }

        /// <summary>
        /// Scales an array. Can be used to scale a vector and a matrix.
        /// </summary>
        /// <param name="alpha">The scalar.</param>
        /// <param name="x">The values to scale.</param>
        /// <param name="result">This result of the scaling.</param>
        /// <remarks>This is similar to the SCAL BLAS routine.</remarks>
        public virtual void ScaleArray(int alpha, int[] x, int[] result)
        {
            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (alpha == 0)
            {
                Array.Clear(result, 0, result.Length);
            }
            else if (alpha == 1)
            {
                x.Copy(result);
            }
            else
            {
                CommonParallel.For(0, x.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        result[i] = alpha * x[i];
                    }
                });
            }
        }

        /// <summary>
        /// Conjugates an array. Can be used to conjugate a vector and a matrix.
        /// </summary>
        /// <param name="x">The values to conjugate.</param>
        /// <param name="result">This result of the conjugation.</param>
        public virtual void ConjugateArray(int[] x, int[] result)
        {
            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (!ReferenceEquals(x, result))
            {
                x.CopyTo(result, 0);
            }
        }

        /// <summary>
        /// Computes the dot product of x and y.
        /// </summary>
        /// <param name="x">The vector x.</param>
        /// <param name="y">The vector y.</param>
        /// <returns>The dot product of x and y.</returns>
        /// <remarks>This is equivalent to the DOT BLAS routine.</remarks>
        public virtual int DotProduct(int[] x, int[] y)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (y.Length != x.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            var sum = 0;

            for (var index = 0; index < y.Length; index++)
            {
                sum += y[index]*x[index];
            }

            return sum;
        }

        /// <summary>
        /// Does a point wise add of two arrays <c>z = x + y</c>. This can be used
        /// to add vectors or matrices.
        /// </summary>
        /// <param name="x">The array x.</param>
        /// <param name="y">The array y.</param>
        /// <param name="result">The result of the addition.</param>
        /// <remarks>There is no equivalent BLAS routine, but many libraries
        /// provide optimized (parallel and/or vectorized) versions of this
        /// routine.</remarks>
        public virtual void AddArrays(int[] x, int[] y, int[] result)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            if (y.Length != x.Length || y.Length != result.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            CommonParallel.For(0, y.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    result[i] = x[i] + y[i];
                }
            });
        }

        /// <summary>
        /// Does a point wise subtraction of two arrays <c>z = x - y</c>. This can be used
        /// to subtract vectors or matrices.
        /// </summary>
        /// <param name="x">The array x.</param>
        /// <param name="y">The array y.</param>
        /// <param name="result">The result of the subtraction.</param>
        /// <remarks>There is no equivalent BLAS routine, but many libraries
        /// provide optimized (parallel and/or vectorized) versions of this
        /// routine.</remarks>
        public virtual void SubtractArrays(int[] x, int[] y, int[] result)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            if (y.Length != x.Length || y.Length != result.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            CommonParallel.For(0, y.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    result[i] = x[i] - y[i];
                }
            });
        }

        /// <summary>
        /// Does a point wise multiplication of two arrays <c>z = x * y</c>. This can be used
        /// to multiple elements of vectors or matrices.
        /// </summary>
        /// <param name="x">The array x.</param>
        /// <param name="y">The array y.</param>
        /// <param name="result">The result of the point wise multiplication.</param>
        /// <remarks>There is no equivalent BLAS routine, but many libraries
        /// provide optimized (parallel and/or vectorized) versions of this
        /// routine.</remarks>
        public virtual void PointWiseMultiplyArrays(int[] x, int[] y, int[] result)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            if (y.Length != x.Length || y.Length != result.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            CommonParallel.For(0, y.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    result[i] = x[i] * y[i];
                }
            });
        }

        /// <summary>
        /// Does a point wise division of two arrays <c>z = x / y</c>. This can be used
        /// to divide elements of vectors or matrices.
        /// </summary>
        /// <param name="x">The array x.</param>
        /// <param name="y">The array y.</param>
        /// <param name="result">The result of the point wise division.</param>
        /// <remarks>There is no equivalent BLAS routine, but many libraries
        /// provide optimized (parallel and/or vectorized) versions of this
        /// routine.</remarks>
        public virtual void PointWiseDivideArrays(int[] x, int[] y, int[] result)
        {
            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            if (y.Length != x.Length || y.Length != result.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            CommonParallel.For(0, y.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    result[i] = x[i] / y[i];
                }
            });
        }

        /// <summary>
        /// Computes the requested <see cref="Norm"/> of the matrix.
        /// </summary>
        /// <param name="norm">The type of norm to compute.</param>
        /// <param name="rows">The number of rows.</param>
        /// <param name="columns">The number of columns.</param>
        /// <param name="matrix">The matrix to compute the norm from.</param>
        /// <returns>
        /// The requested <see cref="Norm"/> of the matrix.
        /// </returns>
        public virtual double MatrixNorm(Norm norm, int rows, int columns, int[] matrix)
        {
            switch (norm)
            {
                case Norm.OneNorm:
                    var norm1 = 0d;
                    for (var j = 0; j < columns; j++)
                    {
                        var s = 0d;
                        for (var i = 0; i < rows; i++)
                        {
                            s += Math.Abs(matrix[(j*rows) + i]);
                        }
                        norm1 = Math.Max(norm1, s);
                    }
                    return norm1;
                case Norm.LargestAbsoluteValue:
                    var normMax = 0d;
                    for (var j = 0; j < columns; j++)
                    {
                        for (var i = 0; i < rows; i++)
                        {
                            normMax = Math.Max(Math.Abs(matrix[(j * rows) + i]), normMax);
                        }
                    }
                    return normMax;
                case Norm.InfinityNorm:
                    var r = new double[rows];
                    for (var j = 0; j < columns; j++)
                    {
                        for (var i = 0; i < rows; i++)
                        {
                            r[i] += Math.Abs(matrix[(j * rows) + i]);
                        }
                    }
                    // TODO: reuse
                    var max = r[0];
                    for (int i = 0; i < r.Length; i++)
                    {
                        if (r[i] > max)
                        {
                            max = r[i];
                        }
                    }
                    return max;
                case Norm.FrobeniusNorm:
                    var aat = new int[rows*rows];
                    MatrixMultiplyWithUpdate(Transpose.DontTranspose, Transpose.Transpose, 1, matrix, rows, columns, matrix, rows, columns, 0, aat);
                    var normF = 0d;
                    for (var i = 0; i < rows; i++)
                    {
                        normF += Math.Abs(aat[(i * rows) + i]);
                    }
                    return Math.Sqrt(normF);
                default:
                    throw new NotSupportedException();
            }
        }

        /// <summary>
        /// Multiples two matrices. <c>result = x * y</c>
        /// </summary>
        /// <param name="x">The x matrix.</param>
        /// <param name="rowsX">The number of rows in the x matrix.</param>
        /// <param name="columnsX">The number of columns in the x matrix.</param>
        /// <param name="y">The y matrix.</param>
        /// <param name="rowsY">The number of rows in the y matrix.</param>
        /// <param name="columnsY">The number of columns in the y matrix.</param>
        /// <param name="result">Where to store the result of the multiplication.</param>
        /// <remarks>This is a simplified version of the BLAS GEMM routine with alpha
        /// set to 1.0 and beta set to 0.0, and x and y are not transposed.</remarks>
        public virtual void MatrixMultiply(int[] x, int rowsX, int columnsX, int[] y, int rowsY, int columnsY, int[] result)
        {
            // First check some basic requirement on the parameters of the matrix multiplication.
            if (x == null)
            {
                throw new ArgumentNullException("x");
            }

            if (y == null)
            {
                throw new ArgumentNullException("y");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            if (rowsX*columnsX != x.Length)
            {
                throw new ArgumentException("x.Length != xRows * xColumns");
            }

            if (rowsY*columnsY != y.Length)
            {
                throw new ArgumentException("y.Length != yRows * yColumns");
            }

            if (columnsX != rowsY)
            {
                throw new ArgumentException("xColumns != yRows");
            }

            if (rowsX*columnsY != result.Length)
            {
                throw new ArgumentException("xRows * yColumns != result.Length");
            }

            // Check whether we will be overwriting any of our inputs and make copies if necessary.
            // TODO - we can don't have to allocate a completely new matrix when x or y point to the same memory
            // as result, we can do it on a row wise basis. We should investigate this.
            int[] xdata;
            if (ReferenceEquals(x, result))
            {
                xdata = (int[]) x.Clone();
            }
            else
            {
                xdata = x;
            }

            int[] ydata;
            if (ReferenceEquals(y, result))
            {
                ydata = (int[]) y.Clone();
            }
            else
            {
                ydata = y;
            }

            Array.Clear(result, 0, result.Length);

            CacheObliviousMatrixMultiply(Transpose.DontTranspose, Transpose.DontTranspose, 1, xdata, 0, 0, ydata, 0, 0, result, 0, 0, rowsX, columnsY, columnsX, rowsX, columnsY, columnsX, true);
        }

        /// <summary>
        /// Multiplies two matrices and updates another with the result. <c>c = alpha*op(a)*op(b) + beta*c</c>
        /// </summary>
        /// <param name="transposeA">How to transpose the <paramref name="a"/> matrix.</param>
        /// <param name="transposeB">How to transpose the <paramref name="b"/> matrix.</param>
        /// <param name="alpha">The value to scale <paramref name="a"/> matrix.</param>
        /// <param name="a">The a matrix.</param>
        /// <param name="rowsA">The number of rows in the <paramref name="a"/> matrix.</param>
        /// <param name="columnsA">The number of columns in the <paramref name="a"/> matrix.</param>
        /// <param name="b">The b matrix</param>
        /// <param name="rowsB">The number of rows in the <paramref name="b"/> matrix.</param>
        /// <param name="columnsB">The number of columns in the <paramref name="b"/> matrix.</param>
        /// <param name="beta">The value to scale the <paramref name="c"/> matrix.</param>
        /// <param name="c">The c matrix.</param>
        public virtual void MatrixMultiplyWithUpdate(Transpose transposeA, Transpose transposeB, int alpha, int[] a, int rowsA, int columnsA, int[] b, int rowsB, int columnsB, int beta, int[] c)
        {
            int m; // The number of rows of matrix op(A) and of the matrix C.
            int n; // The number of columns of matrix op(B) and of the matrix C.
            int k; // The number of columns of matrix op(A) and the rows of the matrix op(B).

            // First check some basic requirement on the parameters of the matrix multiplication.
            if (a == null)
            {
                throw new ArgumentNullException("a");
            }

            if (b == null)
            {
                throw new ArgumentNullException("b");
            }

            if (transposeA > Transpose.DontTranspose && transposeB > Transpose.DontTranspose)
            {
                if (rowsA != columnsB)
                {
                    throw new ArgumentOutOfRangeException();
                }

                if (columnsA*rowsB != c.Length)
                {
                    throw new ArgumentOutOfRangeException();
                }

                m = columnsA;
                n = rowsB;
                k = rowsA;
            }
            else if (transposeA > Transpose.DontTranspose)
            {
                if (rowsA != rowsB)
                {
                    throw new ArgumentOutOfRangeException();
                }

                if (columnsA*columnsB != c.Length)
                {
                    throw new ArgumentOutOfRangeException();
                }

                m = columnsA;
                n = columnsB;
                k = rowsA;
            }
            else if (transposeB > Transpose.DontTranspose)
            {
                if (columnsA != columnsB)
                {
                    throw new ArgumentOutOfRangeException();
                }

                if (rowsA*rowsB != c.Length)
                {
                    throw new ArgumentOutOfRangeException();
                }

                m = rowsA;
                n = rowsB;
                k = columnsA;
            }
            else
            {
                if (columnsA != rowsB)
                {
                    throw new ArgumentOutOfRangeException();
                }

                if (rowsA*columnsB != c.Length)
                {
                    throw new ArgumentOutOfRangeException();
                }

                m = rowsA;
                n = columnsB;
                k = columnsA;
            }

            if (alpha == 0 && beta == 0)
            {
                Array.Clear(c, 0, c.Length);
                return;
            }

            // Check whether we will be overwriting any of our inputs and make copies if necessary.
            // TODO - we can don't have to allocate a completely new matrix when x or y point to the same memory
            // as result, we can do it on a row wise basis. We should investigate this.
            int[] adata;
            if (ReferenceEquals(a, c))
            {
                adata = (int[]) a.Clone();
            }
            else
            {
                adata = a;
            }

            int[] bdata;
            if (ReferenceEquals(b, c))
            {
                bdata = (int[]) b.Clone();
            }
            else
            {
                bdata = b;
            }

            if (beta == 0)
            {
                Array.Clear(c, 0, c.Length);
            }
            else if (beta != 1)
            {
                ScaleArray(beta, c, c);
            }

            if (alpha == 0)
            {
                return;
            }

            CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, adata, 0, 0, bdata, 0, 0, c, 0, 0, m, n, k, m, n, k, true);
        }

        /// <summary>
        /// Cache-Oblivious Matrix Multiplication
        /// </summary>
        /// <param name="transposeA">if set to <c>true</c> transpose matrix A.</param>
        /// <param name="transposeB">if set to <c>true</c> transpose matrix B.</param>
        /// <param name="alpha">The value to scale the matrix A with.</param>
        /// <param name="matrixA">The matrix A.</param>
        /// <param name="shiftArow">Row-shift of the left matrix</param>
        /// <param name="shiftAcol">Column-shift of the left matrix</param>
        /// <param name="matrixB">The matrix B.</param>
        /// <param name="shiftBrow">Row-shift of the right matrix</param>
        /// <param name="shiftBcol">Column-shift of the right matrix</param>
        /// <param name="result">The matrix C.</param>
        /// <param name="shiftCrow">Row-shift of the result matrix</param>
        /// <param name="shiftCcol">Column-shift of the result matrix</param>
        /// <param name="m">The number of rows of matrix op(A) and of the matrix C.</param>
        /// <param name="n">The number of columns of matrix op(B) and of the matrix C.</param>
        /// <param name="k">The number of columns of matrix op(A) and the rows of the matrix op(B).</param>
        /// <param name="constM">The constant number of rows of matrix op(A) and of the matrix C.</param>
        /// <param name="constN">The constant number of columns of matrix op(B) and of the matrix C.</param>
        /// <param name="constK">The constant number of columns of matrix op(A) and the rows of the matrix op(B).</param>
        /// <param name="first">Indicates if this is the first recursion.</param>
        static void CacheObliviousMatrixMultiply(Transpose transposeA, Transpose transposeB, int alpha, int[] matrixA, int shiftArow, int shiftAcol, int[] matrixB, int shiftBrow, int shiftBcol, int[] result, int shiftCrow, int shiftCcol, int m, int n, int k, int constM, int constN, int constK, bool first)
        {
            if (m + n <= Control.ParallelizeOrder)
            {
                if (transposeA > Transpose.DontTranspose && transposeB > Transpose.DontTranspose)
                {
                    for (var m1 = 0; m1 < m; m1++)
                    {
                        var matArowPos = m1 + shiftArow;
                        var matCrowPos = m1 + shiftCrow;
                        for (var n1 = 0; n1 < n; ++n1)
                        {
                            var matBcolPos = n1 + shiftBcol;
                            int sum = 0;
                            for (var k1 = 0; k1 < k; ++k1)
                            {
                                sum += matrixA[(matArowPos*constK) + k1 + shiftAcol]*
                                    matrixB[((k1 + shiftBrow)*constN) + matBcolPos];
                            }

                            result[((n1 + shiftCcol)*constM) + matCrowPos] += alpha*sum;
                        }
                    }
                }
                else if (transposeA > Transpose.DontTranspose)
                {
                    for (var m1 = 0; m1 < m; m1++)
                    {
                        var matArowPos = m1 + shiftArow;
                        var matCrowPos = m1 + shiftCrow;
                        for (var n1 = 0; n1 < n; ++n1)
                        {
                            var matBcolPos = n1 + shiftBcol;
                            int sum = 0;
                            for (var k1 = 0; k1 < k; ++k1)
                            {
                                sum += matrixA[(matArowPos*constK) + k1 + shiftAcol]*
                                    matrixB[(matBcolPos*constK) + k1 + shiftBrow];
                            }

                            result[((n1 + shiftCcol)*constM) + matCrowPos] += alpha*sum;
                        }
                    }
                }
                else if (transposeB > Transpose.DontTranspose)
                {
                    for (var m1 = 0; m1 < m; m1++)
                    {
                        var matArowPos = m1 + shiftArow;
                        var matCrowPos = m1 + shiftCrow;
                        for (var n1 = 0; n1 < n; ++n1)
                        {
                            var matBcolPos = n1 + shiftBcol;
                            int sum = 0;
                            for (var k1 = 0; k1 < k; ++k1)
                            {
                                sum += matrixA[((k1 + shiftAcol)*constM) + matArowPos]*
                                    matrixB[((k1 + shiftBrow)*constN) + matBcolPos];
                            }

                            result[((n1 + shiftCcol)*constM) + matCrowPos] += alpha*sum;
                        }
                    }
                }
                else
                {
                    for (var m1 = 0; m1 < m; m1++)
                    {
                        var matArowPos = m1 + shiftArow;
                        var matCrowPos = m1 + shiftCrow;
                        for (var n1 = 0; n1 < n; ++n1)
                        {
                            var matBcolPos = n1 + shiftBcol;
                            int sum = 0;
                            for (var k1 = 0; k1 < k; ++k1)
                            {
                                sum += matrixA[((k1 + shiftAcol)*constM) + matArowPos]*
                                    matrixB[(matBcolPos*constK) + k1 + shiftBrow];
                            }

                            result[((n1 + shiftCcol)*constM) + matCrowPos] += alpha*sum;
                        }
                    }
                }
            }
            else
            {
                // divide and conquer
                int m2 = m/2, n2 = n/2, k2 = k/2;

                if (first)
                {
                    CommonParallel.Invoke(
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol, matrixB, shiftBrow, shiftBcol, result, shiftCrow, shiftCcol, m2, n2, k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol, matrixB, shiftBrow, shiftBcol + n2, result, shiftCrow, shiftCcol + n2, m2, n - n2, k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol, matrixB, shiftBrow, shiftBcol, result, shiftCrow + m2, shiftCcol, m - m2, n2, k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol, matrixB, shiftBrow, shiftBcol + n2, result, shiftCrow + m2, shiftCcol + n2, m - m2, n - n2, k2, constM, constN, constK, false));

                    CommonParallel.Invoke(
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol, result, shiftCrow, shiftCcol, m2, n2, k - k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol + n2, result, shiftCrow, shiftCcol + n2, m2, n - n2, k - k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol, result, shiftCrow + m2, shiftCcol, m - m2, n2, k - k2, constM, constN, constK, false),
                        () => CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol + n2, result, shiftCrow + m2, shiftCcol + n2, m - m2, n - n2, k - k2, constM, constN, constK, false));
                }
                else
                {
                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol, matrixB, shiftBrow, shiftBcol, result, shiftCrow, shiftCcol, m2, n2, k2, constM, constN, constK, false);
                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol, matrixB, shiftBrow, shiftBcol + n2, result, shiftCrow, shiftCcol + n2, m2, n - n2, k2, constM, constN, constK, false);

                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol, result, shiftCrow, shiftCcol, m2, n2, k - k2, constM, constN, constK, false);
                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol + n2, result, shiftCrow, shiftCcol + n2, m2, n - n2, k - k2, constM, constN, constK, false);

                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol, matrixB, shiftBrow, shiftBcol, result, shiftCrow + m2, shiftCcol, m - m2, n2, k2, constM, constN, constK, false);
                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol, matrixB, shiftBrow, shiftBcol + n2, result, shiftCrow + m2, shiftCcol + n2, m - m2, n - n2, k2, constM, constN, constK, false);

                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol, result, shiftCrow + m2, shiftCcol, m - m2, n2, k - k2, constM, constN, constK, false);
                    CacheObliviousMatrixMultiply(transposeA, transposeB, alpha, matrixA, shiftArow + m2, shiftAcol + k2, matrixB, shiftBrow + k2, shiftBcol + n2, result, shiftCrow + m2, shiftCcol + n2, m - m2, n - n2, k - k2, constM, constN, constK, false);
                }
            }
        }

        /// <summary>
        /// Computes the LUP factorization of A. P*A = L*U.
        /// </summary>
        /// <param name="data">An <paramref name="order"/> by <paramref name="order"/> matrix. The matrix is overwritten with the
        /// the LU factorization on exit. The lower triangular factor L is stored in under the diagonal of <paramref name="data"/> (the diagonal is always 1.0
        /// for the L factor). The upper triangular factor U is stored on and above the diagonal of <paramref name="data"/>.</param>
        /// <param name="order">The order of the square matrix <paramref name="data"/>.</param>
        /// <param name="ipiv">On exit, it contains the pivot indices. The size of the array must be <paramref name="order"/>.</param>
        /// <remarks>This is equivalent to the GETRF LAPACK routine.</remarks>
        public virtual void LUFactor(int[] data, int order, int[] ipiv)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the inverse of matrix using LU factorization.
        /// </summary>
        /// <param name="a">The N by N matrix to invert. Contains the inverse On exit.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <remarks>This is equivalent to the GETRF and GETRI LAPACK routines.</remarks>
        public virtual void LUInverse(int[] a, int order)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the inverse of a previously factored matrix.
        /// </summary>
        /// <param name="a">The LU factored N by N matrix.  Contains the inverse On exit.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <param name="ipiv">The pivot indices of <paramref name="a"/>.</param>
        /// <remarks>This is equivalent to the GETRI LAPACK routine.</remarks>
        public virtual void LUInverseFactored(int[] a, int order, int[] ipiv)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the inverse of matrix using LU factorization.
        /// </summary>
        /// <param name="a">The N by N matrix to invert. Contains the inverse On exit.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <param name="work">The work array. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <remarks>This is equivalent to the GETRF and GETRI LAPACK routines.</remarks>
        public virtual void LUInverse(int[] a, int order, int[] work)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the inverse of a previously factored matrix.
        /// </summary>
        /// <param name="a">The LU factored N by N matrix.  Contains the inverse On exit.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <param name="ipiv">The pivot indices of <paramref name="a"/>.</param>
        /// <param name="work">The work array. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <remarks>This is equivalent to the GETRI LAPACK routine.</remarks>
        public virtual void LUInverseFactored(int[] a, int order, int[] ipiv, int[] work)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using LU factorization.
        /// </summary>
        /// <param name="columnsOfB">The number of columns of B.</param>
        /// <param name="a">The square matrix A.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <param name="b">On entry the B matrix; on exit the X matrix.</param>
        /// <remarks>This is equivalent to the GETRF and GETRS LAPACK routines.</remarks>
        public virtual void LUSolve(int columnsOfB, int[] a, int order, int[] b)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using a previously factored A matrix.
        /// </summary>
        /// <param name="columnsOfB">The number of columns of B.</param>
        /// <param name="a">The factored A matrix.</param>
        /// <param name="order">The order of the square matrix <paramref name="a"/>.</param>
        /// <param name="ipiv">The pivot indices of <paramref name="a"/>.</param>
        /// <param name="b">On entry the B matrix; on exit the X matrix.</param>
        /// <remarks>This is equivalent to the GETRS LAPACK routine.</remarks>
        public virtual void LUSolveFactored(int columnsOfB, int[] a, int order, int[] ipiv, int[] b)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the Cholesky factorization of A.
        /// </summary>
        /// <param name="a">On entry, a square, positive definite matrix. On exit, the matrix is overwritten with the
        /// the Cholesky factorization.</param>
        /// <param name="order">The number of rows or columns in the matrix.</param>
        /// <remarks>This is equivalent to the POTRF LAPACK routine.</remarks>
        public virtual void CholeskyFactor(int[] a, int order)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using Cholesky factorization.
        /// </summary>
        /// <param name="a">The square, positive definite matrix A.</param>
        /// <param name="orderA">The number of rows and columns in A.</param>
        /// <param name="b">On entry the B matrix; on exit the X matrix.</param>
        /// <param name="columnsB">The number of columns in the B matrix.</param>
        /// <remarks>This is equivalent to the POTRF add POTRS LAPACK routines.</remarks>
        public virtual void CholeskySolve(int[] a, int orderA, int[] b, int columnsB)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using a previously factored A matrix.
        /// </summary>
        /// <param name="a">The square, positive definite matrix A. Has to be different than <paramref name="b"/>.</param>
        /// <param name="orderA">The number of rows and columns in A.</param>
        /// <param name="b">On entry the B matrix; on exit the X matrix.</param>
        /// <param name="columnsB">The number of columns in the B matrix.</param>
        /// <remarks>This is equivalent to the POTRS LAPACK routine.</remarks>
        public virtual void CholeskySolveFactored(int[] a, int orderA, int[] b, int columnsB)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the QR factorization of A.
        /// </summary>
        /// <param name="r">On entry, it is the M by N A matrix to factor. On exit,
        /// it is overwritten with the R matrix of the QR factorization. </param>
        /// <param name="rowsR">The number of rows in the A matrix.</param>
        /// <param name="columnsR">The number of columns in the A matrix.</param>
        /// <param name="q">On exit, A M by M matrix that holds the Q matrix of the
        /// QR factorization.</param>
        /// <param name="tau">A min(m,n) vector. On exit, contains additional information
        /// to be used by the QR solve routine.</param>
        /// <remarks>This is similar to the GEQRF and ORGQR LAPACK routines.</remarks>
        public virtual void QRFactor(int[] r, int rowsR, int columnsR, int[] q, int[] tau)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the QR factorization of A.
        /// </summary>
        /// <param name="r">On entry, it is the M by N A matrix to factor. On exit,
        /// it is overwritten with the R matrix of the QR factorization. </param>
        /// <param name="rowsR">The number of rows in the A matrix.</param>
        /// <param name="columnsR">The number of columns in the A matrix.</param>
        /// <param name="q">On exit, A M by M matrix that holds the Q matrix of the
        /// QR factorization.</param>
        /// <param name="tau">A min(m,n) vector. On exit, contains additional information
        /// to be used by the QR solve routine.</param>
        /// <param name="work">The work array. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <remarks>This is similar to the GEQRF and ORGQR LAPACK routines.</remarks>
        public virtual void QRFactor(int[] r, int rowsR, int columnsR, int[] q, int[] tau, int[] work)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the QR factorization of A.
        /// </summary>
        /// <param name="a">On entry, it is the M by N A matrix to factor. On exit,
        /// it is overwritten with the Q matrix of the QR factorization.</param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="r">On exit, A N by N matrix that holds the R matrix of the
        /// QR factorization.</param>
        /// <param name="tau">A min(m,n) vector. On exit, contains additional information
        /// to be used by the QR solve routine.</param>
        /// <remarks>This is similar to the GEQRF and ORGQR LAPACK routines.</remarks>
        public virtual void ThinQRFactor(int[] a, int rowsA, int columnsA, int[] r, int[] tau)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the thin QR factorization of A where M &gt; N.
        /// </summary>
        /// <param name="a">On entry, it is the M by N A matrix to factor. On exit,
        /// it is overwritten with the Q matrix of the QR factorization.</param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="r">On exit, A N by N matrix that holds the R matrix of the
        /// QR factorization.</param>
        /// <param name="tau">A min(m,n) vector. On exit, contains additional information
        /// to be used by the QR solve routine.</param>
        /// <param name="work">The work array. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <remarks>This is similar to the GEQRF and ORGQR LAPACK routines.</remarks>
        public virtual void ThinQRFactor(int[] a, int rowsA, int columnsA, int[] r, int[] tau, int[] work)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using QR factorization of A.
        /// </summary>
        /// <param name="a">The A matrix.</param>
        /// <param name="rows">The number of rows in the A matrix.</param>
        /// <param name="columns">The number of columns in the A matrix.</param>
        /// <param name="b">The B matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        /// <param name="method">The type of QR factorization to perform. <seealso cref="QRMethod"/></param>
        /// <remarks>Rows must be greater or equal to columns.</remarks>
        public virtual void QRSolve(int[] a, int rows, int columns, int[] b, int columnsB, int[] x, QRMethod method = QRMethod.Full)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using QR factorization of A.
        /// </summary>
        /// <param name="a">The A matrix.</param>
        /// <param name="rows">The number of rows in the A matrix.</param>
        /// <param name="columns">The number of columns in the A matrix.</param>
        /// <param name="b">The B matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        /// <param name="work">The work array. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <param name="method">The type of QR factorization to perform. <seealso cref="QRMethod"/></param>
        /// <remarks>Rows must be greater or equal to columns.</remarks>
        public virtual void QRSolve(int[] a, int rows, int columns, int[] b, int columnsB, int[] x, int[] work, QRMethod method = QRMethod.Full)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using a previously QR factored matrix.
        /// </summary>
        /// <param name="q">The Q matrix obtained by QR factor. This is only used for the managed provider and can be
        /// <c>null</c> for the native provider. The native provider uses the Q portion stored in the R matrix.</param>
        /// <param name="r">The R matrix obtained by calling <see cref="QRFactor(int[],int,int,int[],int[])"/>. </param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="tau">Contains additional information on Q. Only used for the native solver
        /// and can be <c>null</c> for the managed provider.</param>
        /// <param name="b">On entry the B matrix; on exit the X matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        /// <param name="work">The work array - only used in the native provider. The array must have a length of at least N,
        /// but should be N*blocksize. The blocksize is machine dependent. On exit, work[0] contains the optimal
        /// work size value.</param>
        /// <param name="method">The type of QR factorization to perform. <seealso cref="QRMethod"/></param>
        /// <remarks>Rows must be greater or equal to columns.</remarks>
        public virtual void QRSolveFactored(int[] q, int[] r, int rowsA, int columnsA, int[] tau, int[] b, int columnsB, int[] x, int[] work, QRMethod method = QRMethod.Full)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using a previously QR factored matrix.
        /// </summary>
        /// <param name="q">The Q matrix obtained by calling <see cref="QRFactor(int[],int,int,int[],int[])"/>.</param>
        /// <param name="r">The R matrix obtained by calling <see cref="QRFactor(int[],int,int,int[],int[])"/>. </param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="tau">Contains additional information on Q. Only used for the native solver
        /// and can be <c>null</c> for the managed provider.</param>
        /// <param name="b">The B matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        /// <param name="method">The type of QR factorization to perform. <seealso cref="QRMethod"/></param>
        /// <remarks>Rows must be greater or equal to columns.</remarks>
        public virtual void QRSolveFactored(int[] q, int[] r, int rowsA, int columnsA, int[] tau, int[] b, int columnsB, int[] x, QRMethod method = QRMethod.Full)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the singular value decomposition of A.
        /// </summary>
        /// <param name="computeVectors">Compute the singular U and VT vectors or not.</param>
        /// <param name="a">On entry, the M by N matrix to decompose. On exit, A may be overwritten.</param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="s">The singular values of A in ascending value.</param>
        /// <param name="u">If <paramref name="computeVectors"/> is <c>true</c>, on exit U contains the left
        /// singular vectors.</param>
        /// <param name="vt">If <paramref name="computeVectors"/> is <c>true</c>, on exit VT contains the transposed
        /// right singular vectors.</param>
        /// <remarks>This is equivalent to the GESVD LAPACK routine.</remarks>
        public virtual void SingularValueDecomposition(bool computeVectors, int[] a, int rowsA, int columnsA, int[] s, int[] u, int[] vt)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the singular value decomposition of A.
        /// </summary>
        /// <param name="computeVectors">Compute the singular U and VT vectors or not.</param>
        /// <param name="a">On entry, the M by N matrix to decompose. On exit, A may be overwritten.</param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="s">The singular values of A in ascending value.</param>
        /// <param name="u">If <paramref name="computeVectors"/> is <c>true</c>, on exit U contains the left
        /// singular vectors.</param>
        /// <param name="vt">If <paramref name="computeVectors"/> is <c>true</c>, on exit VT contains the transposed
        /// right singular vectors.</param>
        /// <param name="work">The work array. Length should be at least <paramref name="rowsA"/>.</param>
        /// <exception cref="NonConvergenceException"></exception>
        public virtual void SingularValueDecomposition(bool computeVectors, int[] a, int rowsA, int columnsA, int[] s, int[] u, int[] vt, int[] work)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using the singular value decomposition of A.
        /// </summary>
        /// <param name="a">On entry, the M by N matrix to decompose.</param>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="b">The B matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        public virtual void SvdSolve(int[] a, int rowsA, int columnsA, int[] b, int columnsB, int[] x)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Solves A*X=B for X using a previously SVD decomposed matrix.
        /// </summary>
        /// <param name="rowsA">The number of rows in the A matrix.</param>
        /// <param name="columnsA">The number of columns in the A matrix.</param>
        /// <param name="s">The s values returned by <see cref="SingularValueDecomposition(bool,int[],int,int,int[],int[],int[])"/>.</param>
        /// <param name="u">The left singular vectors returned by  <see cref="SingularValueDecomposition(bool,int[],int,int,int[],int[],int[])"/>.</param>
        /// <param name="vt">The right singular  vectors returned by  <see cref="SingularValueDecomposition(bool,int[],int,int,int[],int[],int[])"/>.</param>
        /// <param name="b">The B matrix.</param>
        /// <param name="columnsB">The number of columns of B.</param>
        /// <param name="x">On exit, the solution matrix.</param>
        public virtual void SvdSolveFactored(int rowsA, int columnsA, int[] s, int[] u, int[] vt, int[] b, int columnsB, int[] x)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }

        /// <summary>
        /// Computes the eigenvalues and eigenvectors of a matrix.
        /// </summary>
        /// <param name="isSymmetric">Whether the matrix is symmetric or not.</param>
        /// <param name="order">The order of the matrix.</param>
        /// <param name="matrix">The matrix to decompose. The lenth of the array must be order * order.</param>
        /// <param name="matrixEv">On output, the matrix contains the eigen vectors. The lenth of the array must be order * order.</param>
        /// <param name="vectorEv">On output, the eigen values (λ) of matrix in ascending value. The length of the arry must <paramref name="order"/>.</param>
        /// <param name="matrixD">On output, the block diagonal eigenvalue matrix. The lenth of the array must be order * order.</param>
        public virtual void EigenDecomp(bool isSymmetric, int order, int[] matrix, int[] matrixEv, Complex[] vectorEv, int[] matrixD)
        {
            throw new NotSupportedException(Resources.NotSupportedForIntegerMatrices);
        }
    }
}
