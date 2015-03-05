// <copyright file="Milu0.cs" company="Math.NET">
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
using MathNet.Numerics.LinearAlgebra.OneBased.Solvers;
using MathNet.Numerics.LinearAlgebra.OneBased.Storage;
using MathNet.Numerics.Properties;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Double.Solvers
{
    /// <summary>
    /// A simple milu(0) preconditioner.
    /// </summary>
    /// <remarks>
    /// Original Fortran code by Youcef Saad (07 January 2004)
    /// </remarks>
    public sealed class MILU0Preconditioner : IPreconditioner<double>
    {
        // Matrix stored in Modified Sparse Row (MSR) format containing the L and U
        // factors together.

        // The diagonal (stored in alu(1:n) ignore position 0) is inverted. Each i-th row of the matrix
        // contains the i-th row of L (excluding the diagonal entry = 1) followed by
        // the i-th row of U.
        private double[] _alu;

        // The row pointers (stored in jlu(1:n+1) ignore position 0 ) and column indices to off-diagonal elements.
        private int[] _jlu;

        // Pointer to the diagonal elements in MSR storage (for faster LU solving) (ignore position 0).
        private int[] _diag;

        /// <param name="modified">Use modified or standard ILU(0)</param>
        public MILU0Preconditioner(bool modified = true)
        {
            UseModified = modified;
        }

        /// <summary>
        /// Gets or sets a value indicating whether to use modified or standard ILU(0).
        /// </summary>
        public bool UseModified { get; set; }

        /// <summary>
        /// Gets a value indicating whether the preconditioner is initialized.
        /// </summary>
        public bool IsInitialized { get; private set; }

        /// <summary>
        /// Initializes the preconditioner and loads the internal data structures.
        /// </summary>
        /// <param name="matrix">The matrix upon which the preconditioner is based. </param>
        /// <exception cref="ArgumentNullException">If <paramref name="matrix"/> is <see langword="null" />.</exception>
        /// <exception cref="ArgumentException">If <paramref name="matrix"/> is not a square or is not an
        /// instance of SparseCompressedRowMatrixStorage.</exception>
        public void Initialize(Matrix1<double> matrix)
        {
            var csr = matrix.Storage as SparseCompressedRowMatrixStorage<double>;
            if (csr == null)
            {
                throw new ArgumentException(Resources.MatrixMustBeSparse, "matrix");
            }

            // Dimension of matrix
            int n = csr.RowCount;
            if (n != csr.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare, "matrix");
            }

            // Original matrix compressed sparse row storage.
            double[] a = csr.Values;
            int[] ja = csr.ColumnIndices;
            int[] ia = csr.RowPointers;
            int vc = csr.ValueCount;

            _alu = new double[vc + 2];
            _jlu = new int[vc + 2];
            _diag = new int[n + 1];

            int code = Compute(n, a, ja, ia, _alu, _jlu, _diag, UseModified);
            if (code > 0)
            {
                throw new NumericalBreakdownException("Zero pivot encountered on row " + code + " during ILU process");
            }

            IsInitialized = true;
        }

        /// <summary>
        /// Approximates the solution to the matrix equation <b>Ax = b</b>.
        /// </summary>
        /// <param name="input">The right hand side vector b.</param>
        /// <param name="result">The left hand side vector x.</param>
        public void Approximate(Vector1<double> input, Vector1<double> result)
        {
            if (_alu == null)
            {
                throw new ArgumentException(Resources.ArgumentMatrixDoesNotExist);
            }

            int n = _diag.Length - 1;       // bacause _diag has "wasted" position 0

            if ((result.Count != input.Count) || (result.Count != n))
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength);
            }

            // Forward solve.
            for (int i = 1; i <= n; i++)
            {
                var inpi = input.At(i);
                for (int k = _jlu[i]; k < _diag[i]; k++)
                {
                    int jluk = _jlu[k];
                    inpi -= _alu[k] * jluk == i ? inpi : result.At(jluk);
                }
                result.At(i, inpi);
            }

            // Backward solve.
            for (int i = n; i > 0; i--)
            {
                var inpi = input.At(i);
                for (int k = _diag[i]; k <= _jlu[i + 1]; k++)
                {
                    int jluk = _jlu[k];
                    inpi -= _alu[k] * jluk == i ? inpi : result.At(jluk);
                }
                result.At(i, _alu[i] * inpi);
            }
        }

        /// <summary>
        /// MILU0 is a simple milu(0) preconditioner.
        /// </summary>
        /// <param name="n">Order of the matrix.</param>
        /// <param name="a">Matrix values in CSR format (input).</param>
        /// <param name="ja">Column indices (input).</param>
        /// <param name="ia">Row pointers (input).</param>
        /// <param name="alu">Matrix values in MSR format (output).</param>
        /// <param name="jlu">Row pointers and column indices (output).</param>
        /// <param name="ju">Pointer to diagonal elements (output).</param>
        /// <param name="modified">True if the modified/MILU algorithm should be used (recommended)</param>
        /// <returns>Returns 0 on success or k > 0 if a zero pivot was encountered at step k.</returns>
        private int Compute(int n, double[] a, int[] ja, int[] ia, double[] alu, int[] jlu, int[] ju, bool modified)
        {
            var iw = new int[n + 1];
            int i;

            // Set initial pointer value.
            int p = n + 2;
            jlu[1] = p;

            //CONSIDER: since 0 is not a valid index, and .NET guarantees the new iw is already full of 0; this loop is redundant.
            // (A really good optimizer will know this and eliminate this loop...)
            // Initialize work vector.
            for (i = 0; i < n; i++)
            {
                iw[i] = 0;
            }

            // The main loop.
            for (i = 1; i <= n; i++)
            {
                int pold = p;

                // Generating row i of L and U.
                int j;
                for (j = ia[i]; j < ia[i + 1]; j++)
                {
                    // Copy row i of A, JA, IA into row i of ALU, JLU (LU matrix).
                    int jcol = ja[j];

                    if (jcol == i)
                    {
                        alu[i] = a[j];
                        iw[jcol] = i;
                        ju[i] = p;
                    }
                    else
                    {
                        alu[p] = a[j];
                        jlu[p] = ja[j];
                        iw[jcol] = p;
                        p = p + 1;
                    }
                }

                jlu[i + 1] = p;

                double s = 0.0;

                int k;
                for (j = pold; j < ju[i]; j++)
                {
                    int jrow = jlu[j];
                    double tl = alu[j] * alu[jrow];
                    alu[j] = tl;

                    // Perform linear combination.
                    for (k = ju[jrow]; k < jlu[jrow + 1]; k++)
                    {
                        int jw = iw[jlu[k]];
                        if (jw != 0)
                        {
                            alu[jw] = alu[jw] - tl * alu[k];
                        }
                        else
                        {
                            // Accumulate fill-in values.
                            s = s + tl * alu[k];
                        }
                    }
                }

                if (modified)
                {
                    alu[i] = alu[i] - s;
                }

                if (alu[i] == 0.0)
                {
                    return i - 1;       // convert from row index [1..n] to step counter
                }

                // Invert and store diagonal element.
                alu[i] = 1.0 / alu[i];

                // Reset pointers in work array.
                iw[i] = -1;
                for (k = pold; k < p; k++)
                {
                    iw[jlu[k]] = 0;
                }
            }

            return -1;
        }
    }
}
