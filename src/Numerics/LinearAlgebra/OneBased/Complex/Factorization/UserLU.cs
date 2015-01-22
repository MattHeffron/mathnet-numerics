// <copyright file="UserLU.cs" company="Math.NET">
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
using MathNet.Numerics.Properties;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Complex.Factorization
{

#if NOSYSNUMERICS
    using Numerics;
#else
    using System.Numerics;
#endif

    /// <summary>
    /// <para>A class which encapsulates the functionality of an LU factorization.</para>
    /// <para>For a matrix A, the LU factorization is a pair of lower triangular matrix L and
    /// upper triangular matrix U so that A = L*U.</para>
    /// </summary>
    /// <remarks>
    /// The computation of the LU factorization is done at construction time.
    /// </remarks>
    internal sealed class UserLU : LU
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="UserLU"/> class. This object will compute the
        /// LU factorization when the constructor is called and cache it's factorization.
        /// </summary>
        /// <param name="matrix">The matrix to factor.</param>
        /// <exception cref="ArgumentNullException">If <paramref name="matrix"/> is <c>null</c>.</exception>
        /// <exception cref="ArgumentException">If <paramref name="matrix"/> is not a square matrix.</exception>
        public static UserLU Create(Matrix1<Complex> matrix)
        {
            if (matrix == null)
            {
                throw new ArgumentNullException("matrix");
            }

            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            // Create an array for the pivot indices.
            var order = matrix.RowCount;
            var factors = matrix.Clone();
            var pivots1 = new int[order + 1];  // Simplify indexing below, just allocate order+1 elements and "waste" the 0 position

            // Initialize the pivot matrix to the identity permutation.
            for (var i = 1; i <= order; i++)
            {
                pivots1[i] = i;
            }

            var vectorLUcolj = new Complex[order + 1];  // Simplify indexing below, just allocate order+1 elements and "waste" the 0 position
            for (var j = 1; j <= order; j++)
            {
                // Make a copy of the j-th column to localize references.
                for (var i = 1; i <= order; i++)
                {
                    vectorLUcolj[i] = factors.At(i, j);
                }

                // Apply previous transformations.
                for (var i = 1; i <= order; i++)
                {
                    var kmax = Math.Min(i, j);
                    var s = Complex.Zero;
                    for (var k = 1; k <= kmax; k++)
                    {
                        s += factors.At(i, k)*vectorLUcolj[k];
                    }

                    vectorLUcolj[i] -= s;
                    factors.At(i, j, vectorLUcolj[i]);
                }

                // Find pivot and exchange if necessary.
                var p = j;
                for (var i = j + 1; i <= order; i++)
                {
                    if (vectorLUcolj[i].Magnitude > vectorLUcolj[p].Magnitude)
                    {
                        p = i;
                    }
                }

                if (p != j)
                {
                    for (var k = 1; k <= order; k++)
                    {
                        var temp = factors.At(p, k);
                        factors.At(p, k, factors.At(j, k));
                        factors.At(j, k, temp);
                    }

                    pivots1[j] = p;
                }

                // Compute multipliers.
                if (j < order & factors.At(j, j) != 0.0)
                {
                    for (var i = j + 1; i <= order; i++)
                    {
                        factors.At(i, j, (factors.At(i, j)/factors.At(j, j)));
                    }
                }
            }

            // "undo" the wasted 0 position we created above.
            int[] pivots = new int[order];
            for (int i = 0; i < order; i++)
            {
                pivots[i] = pivots1[i + 1] - 1;     // adjust these to be 0-based indices as well.
            }
            return new UserLU(factors, pivots);
        }

        UserLU(Matrix1<Complex> factors, int[] pivots)
            : base(factors, pivots)
        {
        }

        /// <summary>
        /// Solves a system of linear equations, <c>AX = B</c>, with A LU factorized.
        /// </summary>
        /// <param name="input">The right hand side <see cref="Matrix{T}"/>, <c>B</c>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <c>X</c>.</param>
        public override void Solve(Matrix1<Complex> input, Matrix1<Complex> result)
        {
            // Check for proper arguments.
            if (input == null)
            {
                throw new ArgumentNullException("input");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            // Check for proper dimensions.
            if (result.RowCount != input.RowCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameRowDimension);
            }

            if (result.ColumnCount != input.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension);
            }

            if (input.RowCount != Factors.RowCount)
            {
                throw Matrix.DimensionsDontMatch<ArgumentException>(input, Factors);
            }

            // Copy the contents of input to result.
            input.CopyTo(result);
            for (var i = 0; i < Pivots.Length; i++)
            {
                if (Pivots[i] == i)
                {
                    continue;
                }

                var p = Pivots[i] + 1;
                for (var j = 1; j <= result.ColumnCount; j++)
                {
                    var temp = result.At(p, j);
                    result.At(p, j, result.At(i, j));
                    result.At(i, j, temp);
                }
            }

            var order = Factors.RowCount;

            // Solve L*Y = P*B
            for (var k = 1; k <= order; k++)
            {
                for (var i = k + 1; i <= order; i++)
                {
                    var fik = Factors.At(i, k);
                    for (var j = 1; j <= result.ColumnCount; j++)
                    {
                        var temp = result.At(k, j)*fik;
                        result.At(i, j, result.At(i, j) - temp);
                    }
                }
            }

            // Solve U*X = Y;
            for (var k = order; k > 0; k--)
            {
                var fkk = Factors.At(k, k);
                for (var j = 1; j <= result.ColumnCount; j++)
                {
                    result.At(k, j, (result.At(k, j)/fkk));
                }

                for (var i = 1; i <= k; i++)
                {
                    var fik = Factors.At(i, k);
                    for (var j = 1; j <= result.ColumnCount; j++)
                    {
                        var temp = result.At(k, j)*fik;
                        result.At(i, j, result.At(i, j) - temp);
                    }
                }
            }
        }

        /// <summary>
        /// Solves a system of linear equations, <c>Ax = b</c>, with A LU factorized.
        /// </summary>
        /// <param name="input">The right hand side vector, <c>b</c>.</param>
        /// <param name="result">The left hand side <see cref="Matrix{T}"/>, <c>x</c>.</param>
        public override void Solve(Vector1<Complex> input, Vector1<Complex> result)
        {
            // Check for proper arguments.
            if (input == null)
            {
                throw new ArgumentNullException("input");
            }

            if (result == null)
            {
                throw new ArgumentNullException("result");
            }

            // Check for proper dimensions.
            if (input.Count != result.Count)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            if (input.Count != Factors.RowCount)
            {
                throw Matrix.DimensionsDontMatch<ArgumentException>(input, Factors);
            }

            // Copy the contents of input to result.
            input.CopyTo(result);
            for (var i = 0; i < Pivots.Length; i++)
            {
                if (Pivots[i] == i)
                {
                    continue;
                }

                var p = Pivots[i] + 1;
                var temp = result.At(p);
                result.At(p, result.At(i));
                result.At(i, temp);
            }

            var order = Factors.RowCount;

            // Solve L*Y = P*B
            for (var k = 1; k <= order; k++)
            {
                var rk = result.At(k);
                for (var i = k + 1; i <= order; i++)
                {
                    result.At(i, result.At(i) - rk*Factors.At(i, k));
                }
            }

            // Solve U*X = Y;
            for (var k = order; k > 0; k--)
            {
                var rk = result.At(k);
                rk /= Factors.At(k, k);
                result.At(k, rk);
                for (var i = 0; i < k; i++)
                {
                    result.At(i, result.At(i) - rk*Factors.At(i, k));
                }
            }
        }

        /// <summary>
        /// Returns the inverse of this matrix. The inverse is calculated using LU decomposition.
        /// </summary>
        /// <returns>The inverse of this matrix.</returns>
        public override Matrix1<Complex> Inverse()
        {
            var order = Factors.RowCount;
            var inverse = Matrix1<Complex>.Build.SameAs(Factors, order, order);
            for (var i = 1; i <= order; i++)
            {
                inverse.At(i, i, Complex.One);
            }

            return Solve(inverse);
        }
    }
}
