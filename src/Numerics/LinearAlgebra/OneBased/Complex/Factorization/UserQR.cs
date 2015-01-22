// <copyright file="UserQR.cs" company="Math.NET">
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
using System.Linq;
using MathNet.Numerics.LinearAlgebra.OneBased.Factorization;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Threading;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Complex.Factorization
{

#if NOSYSNUMERICS
    using Numerics;
#else
    using System.Numerics;
#endif

    /// <summary>
    /// <para>A class which encapsulates the functionality of the QR decomposition.</para>
    /// <para>Any real square matrix A may be decomposed as A = QR where Q is an orthogonal matrix 
    /// (its columns are orthogonal unit vectors meaning QTQ = I) and R is an upper triangular matrix 
    /// (also called right triangular matrix).</para>
    /// </summary>
    /// <remarks>
    /// The computation of the QR decomposition is done at construction time by Householder transformation.
    /// </remarks>
    internal sealed class UserQR : QR
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="UserQR"/> class. This object will compute the
        /// QR factorization when the constructor is called and cache it's factorization.
        /// </summary>
        /// <param name="matrix">The matrix to factor.</param>
        /// <param name="method">The QR factorization method to use.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="matrix"/> is <c>null</c>.</exception>
        public static UserQR Create(Matrix1<Complex> matrix, QRMethod method = QRMethod.Full)
        {
            if (matrix.RowCount < matrix.ColumnCount)
            {
                throw Matrix.DimensionsDontMatch<ArgumentException>(matrix);
            }

            Matrix1<Complex> q;
            Matrix1<Complex> r;

            ////var minmn = Math.Min(matrix.RowCount, matrix.ColumnCount);
            var minmn = matrix.ColumnCount;         //CONSIDER: the dimension check above, guarantees that ColumnCount == the minimum
            var u = new Complex[minmn + 1][];       // Simplify indexing below, just allocate one extra element and "waste" the 0 position

            if (method == QRMethod.Full)
            {
                r = matrix.Clone();
                int qDim = matrix.RowCount;
                q = Matrix1<Complex>.Build.SameAs(matrix, qDim, qDim);

                for (var i = 1; i <= qDim; i++)
                {
                    q.At(i, i, Complex.One);
                }

                for (var i = 1; i <= minmn; i++)
                {
                    u[i] = GenerateColumn(r, i, i);
                    ComputeQR(u[i], r, i, matrix.RowCount, i + 1, matrix.ColumnCount, Control.MaxDegreeOfParallelism);
                }

                for (var i = minmn; i > 0; i--)
                {
                    ComputeQR(u[i], q, i, qDim, i, qDim, Control.MaxDegreeOfParallelism);
                }
            }
            else
            {
                q = matrix.Clone();

                for (var i = 1; i <= minmn; i++)
                {
                    u[i] = GenerateColumn(q, i, i);
                    ComputeQR(u[i], q, i, q.RowCount, i + 1, q.ColumnCount, Control.MaxDegreeOfParallelism);
                }

                r = q.SubMatrix(1, q.ColumnCount, 1, q.ColumnCount);
                q.Clear();

                for (var i = 1; i <= q.ColumnCount; i++)
                {
                    q.At(i, i, Complex.One);
                }

                for (var i = minmn; i > 0; i--)
                {
                    ComputeQR(u[i], q, i, q.RowCount, i, q.ColumnCount, Control.MaxDegreeOfParallelism);
                }
            }

            return new UserQR(q, r, method);
        }

        UserQR(Matrix1<Complex> q, Matrix1<Complex> rFull, QRMethod method)
            : base(q, rFull, method)
        {
        }

        /// <summary>
        /// Generate column from initial matrix to work array
        /// </summary>
        /// <param name="a">Initial matrix</param>
        /// <param name="row">The first row</param>
        /// <param name="column">Column index</param>
        /// <returns>Generated vector</returns>
        static Complex[] GenerateColumn(Matrix1<Complex> a, int row, int column)
        {
            var ru = a.RowCount - row + 1;      // correct this count (ru) since row parameter is one based
            var u = new Complex[ru];

            for (int i = row; i <= a.RowCount; i++)
            {
                u[i - row] = a.At(i, column);
                a.At(i, column, Complex.Zero);
            }

            var norm = u.Aggregate(Complex.Zero, (current, t) => current + (t.Magnitude*t.Magnitude));
            norm = norm.SquareRoot();

            if (row == a.RowCount || norm.Magnitude == 0)
            {
                a.At(row, column, -u[0]);
                u[0] = Constants.Sqrt2;
                return u;
            }

            if (u[0].Magnitude != 0.0)
            {
                norm = norm.Magnitude*(u[0]/u[0].Magnitude);
            }

            a.At(row, column, -norm);

            for (var i = 0; i < ru; i++)
            {
                u[i] /= norm;
            }

            u[0] += Complex.One;

            var s = u[0].Reciprocal().SquareRoot();
            for (var i = 0; i < ru; i++)
            {
                u[i] = u[i].Conjugate()*s;
            }

            return u;
        }

        /// <summary>
        /// Perform calculation of Q or R
        /// </summary>
        /// <param name="u">Work array</param>
        /// <param name="a">Q or R matrices</param>
        /// <param name="rowStart">The first row</param>
        /// <param name="rowEnd">The last row</param>
        /// <param name="columnStart">The first column</param>
        /// <param name="columnEnd">The last column</param>
        /// <param name="availableCores">Number of available CPUs</param>
        static void ComputeQR(Complex[] u, Matrix1<Complex> a, int rowStart, int rowEnd, int columnStart, int columnEnd, int availableCores)
        {
            if (rowEnd <= rowStart || columnEnd <= columnStart)
            {
                return;
            }

            var tmpColCount = columnEnd - columnStart + 1;

            if ((availableCores > 1) && (tmpColCount > 200))
            {
                var tmpSplit = columnStart + (tmpColCount/2);
                var tmpCores = availableCores/2;

                CommonParallel.Invoke(
                    () => ComputeQR(u, a, rowStart, rowEnd, columnStart, tmpSplit - 1, tmpCores),
                    () => ComputeQR(u, a, rowStart, rowEnd, tmpSplit, columnEnd, tmpCores));
            }
            else
            {
                for (var j = columnStart; j <= columnEnd; j++)
                {
                    var scale = Complex.Zero;
                    for (var i = rowStart; i <= rowEnd; i++)
                    {
                        scale += u[i - rowStart]*a.At(i, j);
                    }

                    for (var i = rowStart; i <= rowEnd; i++)
                    {
                        a.At(i, j, a.At(i, j) - (u[i - rowStart].Conjugate()*scale));
                    }
                }
            }
        }

        /// <summary>
        /// Solves a system of linear equations, <b>AX = B</b>, with A QR factorized.
        /// </summary>
        /// <param name="input">The right hand side <see cref="Matrix{T}"/>, <b>B</b>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <b>X</b>.</param>
        public override void Solve(Matrix1<Complex> input, Matrix1<Complex> result)
        {
            // The solution X should have the same number of columns as B
            if (input.ColumnCount != result.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            }

            // The dimension compatibility conditions for X = A\B require the two matrices A and B to have the same number of rows
            if (FullR.RowCount != input.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameRowDimension);
            }

            // The solution X row dimension is equal to the column dimension of A
            if (FullR.ColumnCount != result.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            }

            var inputCopy = input.Clone();

            // Compute Y = transpose(Q)*B
            var column = new Complex[FullR.RowCount + 1];   // Simplify indexing below, just allocate an extra element and "waste" the 0 position
            for (var j = 1; j <= input.ColumnCount; j++)
            {
                for (var k = 1; k <= FullR.RowCount; k++)
                {
                    column[k] = inputCopy.At(k, j);
                }

                for (var i = 1; i <= FullR.RowCount; i++)
                {
                    var s = Complex.Zero;
                    for (var k = 1; k <= FullR.RowCount; k++)
                    {
                        s += Q.At(k, i).Conjugate()*column[k];
                    }

                    inputCopy.At(i, j, s);
                }
            }

            // Solve R*X = Y;
            for (var k = FullR.ColumnCount; k > 0; k--)
            {
                var frkk = FullR.At(k, k);
                for (var j = 1; j <= input.ColumnCount; j++)
                {
                    inputCopy.At(k, j, inputCopy.At(k, j)/frkk);
                }

                for (var i = 1; i <= k; i++)
                {
                    var frik = FullR.At(i, k);
                    for (var j = 1; j <= input.ColumnCount; j++)
                    {
                        inputCopy.At(i, j, inputCopy.At(i, j) - (inputCopy.At(k, j)*frik));
                    }
                }
            }

            for (var i = 1; i <= FullR.ColumnCount; i++)
            {
                for (var j = 1; j <= inputCopy.ColumnCount; j++)
                {
                    result.At(i, j, inputCopy.At(i, j));
                }
            }
        }

        /// <summary>
        /// Solves a system of linear equations, <b>Ax = b</b>, with A QR factorized.
        /// </summary>
        /// <param name="input">The right hand side vector, <b>b</b>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <b>x</b>.</param>
        public override void Solve(Vector1<Complex> input, Vector1<Complex> result)
        {
            // Ax=b where A is an m x n matrix
            // Check that b is a column vector with m entries
            if (FullR.RowCount != input.Count)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            // Check that x is a column vector with n entries
            if (FullR.ColumnCount != result.Count)
            {
                throw Matrix.DimensionsDontMatch<ArgumentException>(FullR, result);
            }

            var inputCopy = input.Clone();

            // Compute Y = transpose(Q)*B
            var column = new Complex[FullR.RowCount + 1];   // Simplify indexing below, just allocate an extra element and "waste" the 0 position
            for (var k = 1; k <= FullR.RowCount; k++)
            {
                column[k] = inputCopy.At(k);
            }

            for (var i = 1; i <= FullR.RowCount; i++)
            {
                var s = Complex.Zero;
                for (var k = 1; k <= FullR.RowCount; k++)
                {
                    s += Q.At(k, i).Conjugate()*column[k];
                }

                inputCopy.At(i, s);
            }

            // Solve R*X = Y;
            for (var k = FullR.ColumnCount; k > 0; k--)
            {
                var iCk = inputCopy.At(k);
                iCk /= FullR.At(k, k);
                inputCopy.At(k, iCk);
                for (var i = 1; i <= k; i++)
                {
                    inputCopy.At(i, inputCopy.At(i) - iCk*FullR.At(i, k));
                }
            }

            for (var i = 1; i <= FullR.ColumnCount; i++)
            {
                result.At(i, inputCopy.At(i));
            }
        }
    }
}
