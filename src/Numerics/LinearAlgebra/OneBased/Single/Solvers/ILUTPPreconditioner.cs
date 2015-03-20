// <copyright file="Ilutp.cs" company="Math.NET">
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
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra.OneBased.Solvers;
using MathNet.Numerics.Properties;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Single.Solvers
{
    /// <summary>
    /// This class performs an Incomplete LU factorization with drop tolerance
    /// and partial pivoting. The drop tolerance indicates which additional entries
    /// will be dropped from the factorized LU matrices.
    /// </summary>
    /// <remarks>
    /// The ILUTP-Mem algorithm was taken from: <br/>
    /// ILUTP_Mem: a Space-Efficient Incomplete LU Preconditioner
    /// <br/>
    /// Tzu-Yi Chen, Department of Mathematics and Computer Science, <br/>
    /// Pomona College, Claremont CA 91711, USA <br/>
    /// Published in: <br/>
    /// Lecture Notes in Computer Science <br/>
    /// Volume 3046 / 2004 <br/>
    /// pp. 20 - 28 <br/>
    /// Algorithm is described in Section 2, page 22
    /// </remarks>
    public sealed class ILUTPPreconditioner : IPreconditioner<float>
    {
        /// <summary>
        /// The default fill level.
        /// </summary>
        public const double DefaultFillLevel = 200.0;

        /// <summary>
        /// The default drop tolerance.
        /// </summary>
        public const double DefaultDropTolerance = 0.0001;

        /// <summary>
        /// The decomposed upper triangular matrix.
        /// </summary>
        SparseMatrix _upper;

        /// <summary>
        /// The decomposed lower triangular matrix.
        /// </summary>
        SparseMatrix _lower;

        /// <summary>
        /// The array containing the pivot values.
        /// </summary>
        /// <remarks>To simplify indexing into this, just allocate an extra element and "waste" the 0 position.</remarks>
        int[] _pivots;

        /// <summary>
        /// The fill level.
        /// </summary>
        double _fillLevel = DefaultFillLevel;

        /// <summary>
        /// The drop tolerance.
        /// </summary>
        double _dropTolerance = DefaultDropTolerance;

        /// <summary>
        /// The pivot tolerance.
        /// </summary>
        double _pivotTolerance;

        /// <summary>
        /// Initializes a new instance of the <see cref="ILUTPPreconditioner"/> class with the default settings.
        /// </summary>
        public ILUTPPreconditioner()
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ILUTPPreconditioner"/> class with the specified settings.
        /// </summary>
        /// <param name="fillLevel">
        /// The amount of fill that is allowed in the matrix. The value is a fraction of 
        /// the number of non-zero entries in the original matrix. Values should be positive.
        /// </param>
        /// <param name="dropTolerance">
        /// The absolute drop tolerance which indicates below what absolute value an entry 
        /// will be dropped from the matrix. A drop tolerance of 0.0 means that no values
        /// will be dropped. Values should always be positive.
        /// </param>
        /// <param name="pivotTolerance">
        /// The pivot tolerance which indicates at what level pivoting will take place. A
        /// value of 0.0 means that no pivoting will take place.
        /// </param>
        public ILUTPPreconditioner(double fillLevel, double dropTolerance, double pivotTolerance)
        {
            if (fillLevel < 0)
            {
                throw new ArgumentOutOfRangeException("fillLevel");
            }

            if (dropTolerance < 0)
            {
                throw new ArgumentOutOfRangeException("dropTolerance");
            }

            if (pivotTolerance < 0)
            {
                throw new ArgumentOutOfRangeException("pivotTolerance");
            }

            _fillLevel = fillLevel;
            _dropTolerance = dropTolerance;
            _pivotTolerance = pivotTolerance;
        }

        /// <summary>
        /// Gets or sets the amount of fill that is allowed in the matrix. The
        /// value is a fraction of the number of non-zero entries in the original
        /// matrix. The standard value is 200.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Values should always be positive and can be higher than 1.0. A value lower
        /// than 1.0 means that the eventual preconditioner matrix will have fewer
        /// non-zero entries than the original matrix. A value higher than 1.0 means that
        /// the eventual preconditioner can have more non-zero values than the original 
        /// matrix.
        /// </para>
        /// <para>
        /// Note that any changes to the <b>FillLevel</b> after creating the preconditioner
        /// will invalidate the created preconditioner and will require a re-initialization of
        /// the preconditioner.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if a negative value is provided.</exception>
        public double FillLevel
        {
            get { return _fillLevel; }
            set
            {
                if (value < 0)
                {
                    throw new ArgumentOutOfRangeException("value");
                }

                _fillLevel = value;
            }
        }

        /// <summary>
        /// Gets or sets the absolute drop tolerance which indicates below what absolute value
        /// an entry will be dropped from the matrix. The standard value is 0.0001.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The values should always be positive and can be larger than 1.0. A low value will
        /// keep more small numbers in the preconditioner matrix. A high value will remove 
        /// more small numbers from the preconditioner matrix.
        /// </para>
        /// <para>
        /// Note that any changes to the <b>DropTolerance</b> after creating the preconditioner
        /// will invalidate the created preconditioner and will require a re-initialization of
        /// the preconditioner.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if a negative value is provided.</exception>
        public double DropTolerance
        {
            get { return _dropTolerance; }
            set
            {
                if (value < 0)
                {
                    throw new ArgumentOutOfRangeException("value");
                }

                _dropTolerance = value;
            }
        }

        /// <summary>
        /// Gets or sets the pivot tolerance which indicates at what level pivoting will
        /// take place. The standard value is 0.0 which means pivoting will never take place.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The pivot tolerance is used to calculate if pivoting is necessary. Pivoting
        /// will take place if any of the values in a row is bigger than the 
        /// diagonal value of that row divided by the pivot tolerance, i.e. pivoting
        /// will take place if <b>row(i,j) > row(i,i) / PivotTolerance</b> for
        /// any <b>j</b> that is not equal to <b>i</b>.
        /// </para>
        /// <para>
        /// Note that any changes to the <b>PivotTolerance</b> after creating the preconditioner
        /// will invalidate the created preconditioner and will require a re-initialization of
        /// the preconditioner.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown if a negative value is provided.</exception>
        public double PivotTolerance
        {
            get { return _pivotTolerance; }
            set
            {
                if (value < 0)
                {
                    throw new ArgumentOutOfRangeException("value");
                }

                _pivotTolerance = value;
            }
        }

        /// <summary>
        /// Returns the upper triagonal matrix that was created during the LU decomposition.
        /// </summary>
        /// <remarks>
        /// This method is used for debugging purposes only and should normally not be used.
        /// </remarks>
        /// <returns>A new matrix containing the upper triagonal elements.</returns>
        internal Matrix<float> UpperTriangle()
        {
            return _upper.Clone();
        }

        /// <summary>
        /// Returns the lower triagonal matrix that was created during the LU decomposition.
        /// </summary>
        /// <remarks>
        /// This method is used for debugging purposes only and should normally not be used.
        /// </remarks>
        /// <returns>A new matrix containing the lower triagonal elements.</returns>
        internal Matrix<float> LowerTriangle()
        {
            return _lower.Clone();
        }

        /// <summary>
        /// Returns the pivot array. This array is not needed for normal use because
        /// the preconditioner will return the solution vector values in the proper order.
        /// </summary>
        /// <remarks>
        /// This method is used for debugging purposes only and should normally not be used.
        /// </remarks>
        /// <returns>The pivot array.</returns>
        internal int[] Pivots()
        {
            var result = new int[_pivots.Length - 1];
            for (var i = 0; i < result.Length; i++)
            {
                result[i] = _pivots[i + 1];
            }

            return result;
        }

        /// <summary>
        /// Initializes the preconditioner and loads the internal data structures.
        /// </summary>
        /// <param name="matrix">
        /// The <see cref="Matrix"/> upon which this preconditioner is based. Note that the 
        /// method takes a general matrix type. However internally the data is stored 
        /// as a sparse matrix. Therefore it is not recommended to pass a dense matrix.
        /// </param>
        /// <exception cref="ArgumentNullException"> If <paramref name="matrix"/> is <see langword="null" />.</exception>
        /// <exception cref="ArgumentException">If <paramref name="matrix"/> is not a square matrix.</exception>
        public void Initialize(Matrix<float> matrix)
        {
            if (matrix == null)
            {
                throw new ArgumentNullException("matrix");
            }

            if (matrix.RowCount != matrix.ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare, "matrix");
            }

            var sparseMatrix = (matrix is SparseMatrix) ? matrix as SparseMatrix : SparseMatrix.OfMatrix(matrix);

            // The creation of the preconditioner follows the following algorithm.
            // spaceLeft = lfilNnz * nnz(A)
            // for i = 1, .. , n
            // {
            //    w = a(i,*)
            //    for j = 1, .. , i - 1
            //    {
            //        if (w(j) != 0)
            //        {
            //            w(j) = w(j) / a(j,j)
            //            if (w(j) < dropTol)
            //            {
            //                w(j) = 0;
            //            }
            //            if (w(j) != 0)
            //            {
            //                w = w - w(j) * U(j,*)
            //            }
            //        }
            //    }
            //
            //    for j = i, .. ,n
            //    {
            //        if w(j) <= dropTol * ||A(i,*)||
            //        {
            //            w(j) = 0
            //        }
            //    }
            //
            //    spaceRow = spaceLeft / (n - i + 1) // Determine the space for this row
            //    lfil = spaceRow / 2  // space for this row of L
            //    l(i,j) = w(j) for j = 1, .. , i -1 // only the largest lfil elements
            //
            //    lfil = spaceRow - nnz(L(i,:))  // space for this row of U
            //    u(i,j) = w(j) for j = i, .. , n // only the largest lfil - 1 elements
            //    w = 0
            //    
            //    if max(U(i,i + 1: n)) > U(i,i) / pivTol then // pivot if necessary
            //    {
            //        pivot by swapping the max and the diagonal entries
            //        Update L, U
            //        Update P
            //    }
            //    spaceLeft = spaceLeft - nnz(L(i,:)) - nnz(U(i,:))
            // }

            int order = sparseMatrix.RowCount;
            // Create the lower triangular matrix
            _lower = new SparseMatrix(order);

            // Create the upper triangular matrix and copy the values
            _upper = new SparseMatrix(order);

            // Create the pivot array
            _pivots = new int[order + 1];       // Simplify indexing below, just allocate an extra element and "waste" the 0 position
            for (var i = 1; i <= order; i++)
            {
                _pivots[i] = i;
            }

            var workVector = new DenseVector(order);
            var rowVector = new DenseVector(order);
            var indexSorting = new int[order];

            // spaceLeft = lfilNnz * nnz(A)
            var spaceLeft = (int) _fillLevel*sparseMatrix.NonZerosCount;

            // for i = 1, .. , n
            for (var i = 1; i <= order; i++)
            {
                // w = a(i,*)
                sparseMatrix.Row(i, workVector);

                // pivot the row
                PivotRow(workVector);
                var vectorNorm = workVector.InfinityNorm();

                // for j = 1, .. , i - 1)
                for (var j = 1; j < i; j++)
                {
                    // if (w(j) != 0)
                    // {
                    //     w(j) = w(j) / a(j,j)
                    //     if (w(j) < dropTol)
                    //     {
                    //         w(j) = 0;
                    //     }
                    //     if (w(j) != 0)
                    //     {
                    //         w = w - w(j) * U(j,*)
                    //     }
                    var wVj = workVector.At(j);
                    if (wVj != 0.0f)
                    {
                        // Calculate the multiplication factors that go into the L matrix
                        workVector.At(j, wVj /= _upper.At(j, j));
                        if (Math.Abs(wVj) < _dropTolerance)
                        {
                            workVector.At(j, wVj = 0.0f);
                        }

                        // Calculate the addition factor
                        if (wVj != 0.0f)
                        {
                            // vector update all in one go
                            _upper.Row(j, rowVector);

                            // zero out columnVector[k] because we don't need that
                            // one anymore for k = 1 to k = j
                            for (var k = 1; k <= j; k++)
                            {
                                rowVector[k] = 0.0f;
                            }

                            rowVector.Multiply(wVj, rowVector);
                            workVector.Subtract(rowVector, workVector);
                        }
                    }
                }

                // for j = i, .. ,n
                for (var j = i; j <= order; j++)
                {
                    // if w(j) <= dropTol * ||A(i,*)||
                    // {
                    //     w(j) = 0
                    // }
                    if (Math.Abs(workVector[j]) <= _dropTolerance*vectorNorm)
                    {
                        workVector[j] = 0.0f;
                    }
                }

                // spaceRow = spaceLeft / (n - i + 1) // Determine the space for this row
                var spaceRow = spaceLeft / (order - i + 2);       // CONSIDER: in 0-based impl. this should NOT have +1, but for consistency, adjust it here.

                // lfil = spaceRow / 2  // space for this row of L
                var fillLevel = spaceRow/2;
                FindLargestItems(1, i - 1, indexSorting, workVector);

                // l(i,j) = w(j) for j = 1, .. , i -1 // only the largest lfil elements
                var lowerNonZeroCount = 0;
                var count = 0;
                for (var j = 0; j < i - 1; j++)     // i is one based indexing and j is 0 based here (tread carefully!)
                {
                    int iSj = indexSorting[j]; // indexSorting is 0-based but contains 1-based index values
                    if ((count > fillLevel) || (iSj == 0))
                    {
                        break;
                    }

                    _lower.At(i, iSj, workVector.At(iSj));
                    count++;
                    lowerNonZeroCount++;
                }

                FindLargestItems(i + 1, order, indexSorting, workVector);

                // lfil = spaceRow - nnz(L(i,:))  // space for this row of U
                fillLevel = spaceRow - lowerNonZeroCount;

                // u(i,j) = w(j) for j = i + 1, .. , n // only the largest lfil - 1 elements
                var upperNonZeroCount = 0;
                count = 0;
                for (var j = 0; j <= order - i; j++)
                {
                    int iSj = indexSorting[j]; // indexSorting is 0-based but contains 1-based index values
                    if ((count > fillLevel - 1) || (iSj == 0))
                    {
                        break;
                    }

                    _upper.At(i, iSj, workVector.At(iSj));
                    count++;
                    upperNonZeroCount++;
                }

                // Simply copy the diagonal element. Next step is to see if we pivot
                _upper.At(i, i, workVector.At(i));

                // if max(U(i,i + 1: n)) > U(i,i) / pivTol then // pivot if necessary
                // {
                //     pivot by swapping the max and the diagonal entries
                //     Update L, U
                //     Update P
                // }

                // Check if we really need to pivot. If (i+1) >= sparseMatrix.RowCount then
                // we are working on the last row. That means that there is only one number
                // And pivoting is useless. Also the indexSorting array will only contain
                // 0 values.
                if ((i + 1) < order)
                {
                    int iS0 = indexSorting[0];
                    if (Math.Abs(workVector[i]) < _pivotTolerance*Math.Abs(workVector[iS0]))
                    {
                        // swap columns of u (which holds the values of A in the
                        // sections that haven't been partitioned yet.
                        SwapColumns(_upper, i, iS0);

                        // Update P
                        var temp = _pivots[i];
                        _pivots[i] = _pivots[iS0];
                        _pivots[iS0] = temp;
                    }
                }

                // spaceLeft = spaceLeft - nnz(L(i,:)) - nnz(U(i,:))
                spaceLeft -= lowerNonZeroCount + upperNonZeroCount;
            }

            for (var i = 1; i <= order; i++)
            {
                _lower.At(i, i, 1.0f);
            }
        }

        /// <summary>
        /// Pivot elements in the <paramref name="row"/> according to internal pivot array
        /// </summary>
        /// <param name="row">Row <see cref="Vector"/> to pivot in</param>
        void PivotRow(Vector<float> row)
        {
            var knownPivots = new Dictionary<int, int>();

            // pivot the row
            for (var i = 1; i <= row.Count; i++)
            {
                int pi = _pivots[i];
                if ((pi != i) && (!PivotMapFound(knownPivots, i)))
                {
                    // store the pivots in the hashtable
                    knownPivots.Add(pi, i);

                    var t = row.At(i);
                    row.At(i, row.At(pi));
                    row.At(pi, t);
                }
            }
        }

        /// <summary>
        /// Was pivoting already performed
        /// </summary>
        /// <param name="knownPivots">Pivots already done</param>
        /// <param name="currentItem">Current item to pivot</param>
        /// <returns><c>true</c> if performed, otherwise <c>false</c></returns>
        bool PivotMapFound(Dictionary<int, int> knownPivots, int currentItem)
        {
            int pci = _pivots[currentItem];
            int knownItem;
            if (knownPivots.TryGetValue(pci, out knownItem))
            {
                if (knownItem.Equals(currentItem))
                {
                    return true;
                }
            }

            if (knownPivots.TryGetValue(currentItem, out knownItem))
            {
                if (knownItem.Equals(pci))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Swap columns in the <see cref="Matrix"/>
        /// </summary>
        /// <param name="matrix">Source <see cref="Matrix"/>.</param>
        /// <param name="firstColumn">First column index to swap</param>
        /// <param name="secondColumn">Second column index to swap</param>
        static void SwapColumns(Matrix<float> matrix, int firstColumn, int secondColumn)
        {
            for (var i = 1; i <= matrix.RowCount; i++)
            {
                var temp = matrix[i, firstColumn];
                matrix[i, firstColumn] = matrix[i, secondColumn];
                matrix[i, secondColumn] = temp;
            }
        }

        /// <summary>
        /// Sort vector descending, not changing vector but placing sorted indicies to <paramref name="sortedIndices"/>
        /// </summary>
        /// <param name="lowerBound">Start sort form</param>
        /// <param name="upperBound">Sort till upper bound</param>
        /// <param name="sortedIndices">Array with sorted vector indicies (array is 0-based, but contains 1-based index values!)</param>
        /// <param name="values">Source <see cref="Vector"/></param>
        static void FindLargestItems(int lowerBound, int upperBound, int[] sortedIndices, Vector<float> values)
        {
            // Copy the indices for the values into the array
            for (var i = 0; i < upperBound + 1 - lowerBound; i++)
            {
                sortedIndices[i] = lowerBound + i;
            }

            for (var i = upperBound + 1 - lowerBound; i < sortedIndices.Length; i++)
            {
                sortedIndices[i] = 0;
            }

            // Sort the first set of items.
            // Sorting starts at index 0 because the index array
            // starts at zero
            // and ends at index upperBound - lowerBound
            ILUTPElementSorter.SortDoubleIndicesDecreasing(0, upperBound - lowerBound, sortedIndices, values);
        }

        /// <summary>
        /// Approximates the solution to the matrix equation <b>Ax = b</b>.
        /// </summary>
        /// <param name="rhs">The right hand side vector.</param>
        /// <param name="lhs">The left hand side vector. Also known as the result vector.</param>
        public void Approximate(Vector<float> rhs, Vector<float> lhs)
        {
            if (_upper == null)
            {
                throw new ArgumentException(Resources.ArgumentMatrixDoesNotExist);
            }

            if ((lhs.Count != rhs.Count) || (lhs.Count != _upper.RowCount))
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength, "rhs");
            }

            // Solve equation here
            // Pivot(vector, result);
            // Solve L*Y = B(piv,:)
            var rowValues = new DenseVector(_lower.RowCount);
            for (var i = 1; i <= _lower.RowCount; i++)
            {
                _lower.Row(i, rowValues);

                var sum = 0.0f;
                for (var j = 1; j <= i; j++)
                {
                    sum += rowValues.At(j)*lhs.At(j);
                }

                lhs.At(i, rhs.At(i) - sum);
            }

            // Solve U*X = Y;
            for (var i = _upper.RowCount; i > 0; i--)
            {
                _upper.Row(i, rowValues);

                var sum = 0.0f;
                for (var j = _upper.RowCount; j > i; j--)
                {
                    sum += rowValues.At(j)*lhs.At(j);
                }

                lhs.At(i, (1/rowValues.At(i))*(lhs.At(i) - sum));
            }

            // We have a column pivot so we only need to pivot the
            // end result not the incoming right hand side vector
            var temp = lhs.Clone();

            Pivot(temp, lhs);
        }

        /// <summary>
        /// Pivot elements in <see cref="Vector"/> according to internal pivot array
        /// </summary>
        /// <param name="vector">Source <see cref="Vector"/>.</param>
        /// <param name="result">Result <see cref="Vector"/> after pivoting.</param>
        void Pivot(Vector<float> vector, Vector<float> result)
        {
            for (var i = 1; i < _pivots.Length; i++)
            {
                result.At(i, vector.At(_pivots[i]));
            }
        }
    }

    /// <summary>
    /// An element sort algorithm for the <see cref="ILUTPPreconditioner"/> class.
    /// </summary>
    /// <remarks>
    /// This sort algorithm is used to sort the columns in a sparse matrix based on
    /// the value of the element on the diagonal of the matrix.
    /// </remarks>
    internal static class ILUTPElementSorter
    {
        /// <summary>
        /// Sorts the elements of the <paramref name="values"/> vector in decreasing
        /// fashion. The vector itself is not affected.
        /// </summary>
        /// <param name="lowerBound">The starting index.</param>
        /// <param name="upperBound">The stopping index.</param>
        /// <param name="sortedIndices">An array that will contain the sorted indices once the algorithm finishes.</param>
        /// <param name="values">The <see cref="Vector{T}"/> that contains the values that need to be sorted.</param>
        public static void SortDoubleIndicesDecreasing(int lowerBound, int upperBound, int[] sortedIndices, Vector<float> values)
        {
            // Move all the indices that we're interested in to the beginning of the
            // array. Ignore the rest of the indices.
            if (lowerBound > 0)
            {
                for (var i = 0; i < (upperBound - lowerBound + 1); i++)
                {
                    Exchange(sortedIndices, i, i + lowerBound);
                }

                upperBound -= lowerBound;
                lowerBound = 0;
            }

            HeapSortDoublesIndices(lowerBound, upperBound, sortedIndices, values);
        }

        /// <summary>
        /// Sorts the elements of the <paramref name="values"/> vector in decreasing
        /// fashion using heap sort algorithm. The vector itself is not affected.
        /// </summary>
        /// <param name="lowerBound">The starting index.</param>
        /// <param name="upperBound">The stopping index.</param>
        /// <param name="sortedIndices">An array that will contain the sorted indices once the algorithm finishes.</param>
        /// <param name="values">The <see cref="Vector{T}"/> that contains the values that need to be sorted.</param>
        private static void HeapSortDoublesIndices(int lowerBound, int upperBound, int[] sortedIndices, Vector<float> values)
        {
            var start = ((upperBound - lowerBound + 1) / 2) - 1 + lowerBound;
            var end = (upperBound - lowerBound + 1) - 1 + lowerBound;

            BuildDoubleIndexHeap(start, upperBound - lowerBound + 1, sortedIndices, values);

            while (end >= lowerBound)
            {
                Exchange(sortedIndices, end, lowerBound);
                SiftDoubleIndices(sortedIndices, values, lowerBound, end);
                end -= 1;
            }
        }

        /// <summary>
        /// Build heap for double indicies
        /// </summary>
        /// <param name="start">Root position</param>
        /// <param name="count">Length of <paramref name="values"/></param>
        /// <param name="sortedIndices">Indicies of <paramref name="values"/></param>
        /// <param name="values">Target <see cref="Vector{T}"/></param>
        private static void BuildDoubleIndexHeap(int start, int count, int[] sortedIndices, Vector<float> values)
        {
            while (start >= 0)
            {
                SiftDoubleIndices(sortedIndices, values, start, count);
                start--;
            }
        }

        /// <summary>
        /// Sift double indicies
        /// </summary>
        /// <param name="sortedIndices">Indicies of <paramref name="values"/></param>
        /// <param name="values">Target <see cref="Vector{T}"/></param>
        /// <param name="begin">Root position</param>
        /// <param name="count">Length of <paramref name="values"/></param>
        private static void SiftDoubleIndices(int[] sortedIndices, Vector<float> values, int begin, int count)
        {
            var root = begin;

            while (root * 2 < count)
            {
                var child = root * 2;
                if ((child < count - 1) && (values[sortedIndices[child]] > values[sortedIndices[child + 1]]))
                {
                    child += 1;
                }

                if (values[sortedIndices[root]] <= values[sortedIndices[child]])
                {
                    return;
                }

                Exchange(sortedIndices, root, child);
                root = child;
            }
        }

        /// <summary>
        /// Sorts the given integers in a decreasing fashion.
        /// </summary>
        /// <param name="values">The values.</param>
        public static void SortIntegersDecreasing(int[] values)
        {
            HeapSortIntegers(values, values.Length);
        }

        /// <summary>
        /// Sort the given integers in a decreasing fashion using heapsort algorithm 
        /// </summary>
        /// <param name="values">Array of values to sort</param>
        /// <param name="count">Length of <paramref name="values"/></param>
        private static void HeapSortIntegers(int[] values, int count)
        {
            var start = (count / 2) - 1;
            var end = count - 1;

            BuildHeap(values, start, count);

            while (end >= 0)
            {
                Exchange(values, end, 0);
                Sift(values, 0, end);
                end -= 1;
            }
        }

        /// <summary>
        /// Build heap
        /// </summary>
        /// <param name="values">Target values array</param>
        /// <param name="start">Root position</param>
        /// <param name="count">Length of <paramref name="values"/></param>
        private static void BuildHeap(int[] values, int start, int count)
        {
            while (start >= 0)
            {
                Sift(values, start, count);
                start -= 1;
            }
        }

        /// <summary>
        /// Sift values
        /// </summary>
        /// <param name="values">Target value array</param>
        /// <param name="start">Root position</param>
        /// <param name="count">Length of <paramref name="values"/></param>
        private static void Sift(int[] values, int start, int count)
        {
            var root = start;

            while (root * 2 < count)
            {
                var child = root * 2;
                if ((child < count - 1) && (values[child] > values[child + 1]))
                {
                    child += 1;
                }

                if (values[root] > values[child])
                {
                    Exchange(values, root, child);
                    root = child;
                }
                else
                {
                    return;
                }
            }
        }

        /// <summary>
        /// Exchange values in array
        /// </summary>
        /// <param name="values">Target values array</param>
        /// <param name="first">First value to exchange</param>
        /// <param name="second">Second value to exchange</param>
        private static void Exchange(int[] values, int first, int second)
        {
            var t = values[first];
            values[first] = values[second];
            values[second] = t;
        }
    }
}
