// <copyright file="AssertHelpers.cs" company="Math.NET">
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

using System;
using MathNet.Numerics.LinearAlgebra.OneBased;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased
{
#if NOSYSNUMERICS
    using Complex64 = Numerics.Complex;
#else
    using Complex64 = System.Numerics.Complex;
#endif

    /// <summary>
    /// Matrix utility functions to simplify tests.
    /// </summary>
    static public class MatrixHelpers
    {
        /// <summary>
        /// Forces a matrix elements to symmetric. Copies the lower triangle to the upper triangle.
        /// </summary>
        /// <typeparam name="T">The matrix type.</typeparam>
        /// <param name="matrix">The matrix to make symmetric.</param>
        static public void ForceSymmetric<T>(Matrix<T> matrix) where T : struct, IEquatable<T>, IFormattable
        {
            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException("matrix must be square.", "matrix");
            }
            for (var row = 1; row <= matrix.RowCount; row++)
            {
                for (var column = 1; column <= row; column++)
                {
                    matrix.At(column, row, matrix.At(row, column));
                }
            }
        }

        /// <summary>
        /// Forces a matrix elements to hermitian (conjugate symmetric). Copies the conjugate of the values
        /// from the lower triangle to the upper triangle.
        /// </summary>
        /// <param name="matrix">The matrix to make conjugate symmetric.</param>
        static public void ForceHermitian(Matrix<Complex64> matrix)
        {
            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException("matrix must be square.", "matrix");
            }
            for (var row = 1; row <= matrix.RowCount; row++)
            {
                for (var column = 1; column <= row; column++)
                {
                    matrix.At(column, row, matrix.At(row, column).Conjugate());
                }
            }
        }

        /// <summary>
        /// Forces a matrix elements to conjugate symmetric. Copies the conjugate of the values
        /// from the lower triangle to the upper triangle.
        /// </summary>
        /// <param name="matrix">The matrix to make conjugate symmetric.</param>
        public static void ForceHermitian(Matrix<Numerics.Complex32> matrix)
        {
            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException("matrix must be square.", "matrix");
            }
            for (var row = 1; row <= matrix.RowCount; row++)
            {
                for (var column = 1; column <= row; column++)
                {
                    matrix.At(column, row, matrix.At(row, column).Conjugate());
                }
            }
        }

        /// <summary>
        /// Create the matrix for solving the Poisson equation
        /// on a rectangular grid.
        /// </summary>
        /// <typeparam name="T">The matrix type.</typeparam>
        /// <param name="gridSize">The gridsize within the matrix.</param>
        static public Matrix<T> MakePoissonTestMatrix<T>(int gridSize) where T : struct, IEquatable<T>, IFormattable
        {
            // Assemble the matrix. We assume we're solving the Poisson equation
            // on a rectangular gridSize x gridSize grid
            int fullSize = gridSize * gridSize;

            // Create the matrix
            var matrix = Matrix<T>.Build.Sparse(fullSize, fullSize);

            // This is a big HACK to get around the fact that I can't constrain T to be convertable from int
            // But since we know HOW this is being used in testing...
            T minusOne = (T)(object)-1;
            T four = (T)(object)4;

            // The pattern is:
            // 0 .... 0 -1 0 0 0 0 0 0 0 0 -1 4 -1 0 0 0 0 0 0 0 0 -1 0 0 ... 0
            for (var i = 1; i <= matrix.RowCount; i++)
            {
                // Insert the first set of -1's
                if (i > gridSize)
                {
                    matrix[i, i - gridSize + 1] = minusOne;
                }

                // Insert the second set of -1's
                if (i > 1)
                {
                    matrix[i, i - 1] = minusOne;
                }

                // Insert the centerline values
                matrix[i, i] = four;

                // Insert the first trailing set of -1's
                if (i < matrix.RowCount)
                {
                    matrix[i, i + 1] = minusOne;
                }

                // Insert the second trailing set of -1's
                if (i < matrix.RowCount - gridSize + 1)
                {
                    matrix[i, i + gridSize - 1] = minusOne;
                }
            }

            return matrix;
        }
    }
}
