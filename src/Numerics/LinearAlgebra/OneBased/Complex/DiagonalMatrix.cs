﻿// <copyright file="DiagonalMatrix.cs" company="Math.NET">
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
using System.Diagnostics;
using System.Linq;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra.OneBased.Storage;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Threading;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Complex
{

#if NOSYSNUMERICS
    using Numerics;
#else
    using System.Numerics;
#endif

    /// <summary>
    /// A matrix type for diagonal matrices.
    /// </summary>
    /// <remarks>
    /// Diagonal matrices can be non-square matrices but the diagonal always starts
    /// at element 1,1. A diagonal matrix will throw an exception if non diagonal
    /// entries are set. The exception to this is when the off diagonal elements are
    /// 0.0 or NaN; these settings will cause no change to the diagonal matrix.
    /// </remarks>
    [Serializable]
    [DebuggerDisplay("DiagonalMatrix[1] {RowCount}x{ColumnCount}-Complex")]
    public class DiagonalMatrix : Matrix
    {
        /// <summary>
        /// Gets the matrix's data.
        /// </summary>
        /// <value>The matrix's data.</value>
        readonly Complex[] _data;

        /// <summary>
        /// Create a new diagonal matrix straight from an initialized matrix storage instance.
        /// The storage is used directly without copying.
        /// Intended for advanced scenarios where you're working directly with
        /// storage for performance or interop reasons.
        /// </summary>
        public DiagonalMatrix(DiagonalMatrixStorage<Complex> storage)
            : base(storage)
        {
            _data = storage.Data;
        }

        /// <summary>
        /// Create a new square diagonal matrix with the given number of rows and columns.
        /// All cells of the matrix will be initialized to zero.
        /// </summary>
        /// <exception cref="ArgumentException">If the order is less than zero.</exception>
        public DiagonalMatrix(int order)
            : this(new DiagonalMatrixStorage<Complex>(order, order))
         {
         }

        /// <summary>
        /// Create a new diagonal matrix with the given number of rows and columns.
        /// All cells of the matrix will be initialized to zero.
        /// </summary>
        /// <exception cref="ArgumentException">If the row or column count is less than zero.</exception>
        public DiagonalMatrix(int rows, int columns)
            : this(new DiagonalMatrixStorage<Complex>(rows, columns))
        {
        }

        /// <summary>
        /// Create a new diagonal matrix with the given number of rows and columns.
        /// All diagonal cells of the matrix will be initialized to the provided value, all non-diagonal ones to zero.
        /// </summary>
        /// <exception cref="ArgumentException">If the row or column count is less than zero.</exception>
        public DiagonalMatrix(int rows, int columns, Complex diagonalValue)
            : this(rows, columns)
        {
            for (var i = 0; i < _data.Length; i++)
            {
                _data[i] = diagonalValue;
            }
        }

        /// <summary>
        /// Create a new diagonal matrix with the given number of rows and columns directly binding to a raw array.
        /// The array is assumed to contain the diagonal elements only and is used directly without copying.
        /// Very efficient, but changes to the array and the matrix will affect each other.
        /// </summary>
        public DiagonalMatrix(int rows, int columns, Complex[] diagonalStorage)
            : this(new DiagonalMatrixStorage<Complex>(rows, columns, diagonalStorage))
        {
        }

        /// <summary>
        /// Create a new diagonal matrix as a copy of the given other matrix.
        /// This new matrix will be independent from the other matrix.
        /// The matrix to copy from must be diagonal as well.
        /// A new memory block will be allocated for storing the matrix.
        /// </summary>
        public static DiagonalMatrix OfMatrix(Matrix<Complex> matrix)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfMatrix(matrix.Storage));
        }

        /// <summary>
        /// Create a new diagonal matrix as a copy of the given two-dimensional array.
        /// This new matrix will be independent from the provided array.
        /// The array to copy from must be diagonal as well.
        /// A new memory block will be allocated for storing the matrix.
        /// </summary>
        public static DiagonalMatrix OfArray(Complex[,] array)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfArray(array));
        }

        /// <summary>
        /// Create a new diagonal matrix and initialize each diagonal value from the provided indexed enumerable.
        /// Keys must be provided at most once, zero is assumed if a key is omitted.
        /// This new matrix will be independent from the enumerable.
        /// A new memory block will be allocated for storing the matrix.
        /// </summary>
        public static DiagonalMatrix OfIndexedDiagonal(int rows, int columns, IEnumerable<Tuple<int, Complex>> diagonal)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfIndexedEnumerable(rows, columns, diagonal));
        }

        /// <summary>
        /// Create a new diagonal matrix and initialize each diagonal value from the provided enumerable.
        /// This new matrix will be independent from the enumerable.
        /// A new memory block will be allocated for storing the matrix.
        /// </summary>
        public static DiagonalMatrix OfDiagonal(int rows, int columns, IEnumerable<Complex> diagonal)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfEnumerable(rows, columns, diagonal));
        }

        /// <summary>
        /// Create a new diagonal matrix and initialize each diagonal value using the provided init function.
        /// </summary>
        public static DiagonalMatrix Create(int rows, int columns, Func<int, Complex> init)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfInit(rows, columns, init));
        }

        /// <summary>
        /// Create a new square sparse identity matrix where each diagonal value is set to One.
        /// </summary>
        public static DiagonalMatrix CreateIdentity(int order)
        {
            return new DiagonalMatrix(DiagonalMatrixStorage<Complex>.OfValue(order, order, One));
        }

        /// <summary>
        /// Create a new diagonal matrix with diagonal values sampled from the provided random distribution.
        /// </summary>
        public static DiagonalMatrix CreateRandom(int rows, int columns, IContinuousDistribution distribution)
        {
            return new DiagonalMatrix(new DiagonalMatrixStorage<Complex>(rows, columns, Generate.RandomComplex(Math.Min(rows, columns), distribution)));
        }

        /// <summary>
        /// Negate each element of this matrix and place the results into the result matrix.
        /// </summary>
        /// <param name="result">The result of the negation.</param>
        protected override void DoNegate(Matrix<Complex> result)
        {
            var diagResult = result as DiagonalMatrix;
            if (diagResult != null)
            {
                Control.LinearAlgebraProvider.ScaleArray(-1, _data, diagResult._data);
                return;
            }

            result.Clear();
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, -_data[i - 1]);
            }
        }

        /// <summary>
        /// Complex conjugates each element of this matrix and place the results into the result matrix.
        /// </summary>
        /// <param name="result">The result of the conjugation.</param>
        protected override void DoConjugate(Matrix<Complex> result)
        {
            var diagResult = result as DiagonalMatrix;
            if (diagResult != null)
            {
                Control.LinearAlgebraProvider.ConjugateArray(_data, diagResult._data);
                return;
            }

            result.Clear();
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, _data[i - 1].Conjugate());
            }
        }

        /// <summary>
        /// Adds another matrix to this matrix.
        /// </summary>
        /// <param name="other">The matrix to add to this matrix.</param>
        /// <param name="result">The matrix to store the result of the addition.</param>
        /// <exception cref="ArgumentOutOfRangeException">If the two matrices don't have the same dimensions.</exception>
        protected override void DoAdd(Matrix<Complex> other, Matrix<Complex> result)
        {
            // diagonal + diagonal = diagonal
            var diagOther = other as DiagonalMatrix;
            var diagResult = result as DiagonalMatrix;
            if (diagOther != null && diagResult != null)
            {
                Control.LinearAlgebraProvider.AddArrays(_data, diagOther._data, diagResult._data);
                return;
            }

            other.CopyTo(result);
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, result.At(i, i) + _data[i - 1]);
            }
        }

        /// <summary>
        /// Subtracts another matrix from this matrix.
        /// </summary>
        /// <param name="other">The matrix to subtract.</param>
        /// <param name="result">The matrix to store the result of the subtraction.</param>
        /// <exception cref="ArgumentOutOfRangeException">If the two matrices don't have the same dimensions.</exception>
        protected override void DoSubtract(Matrix<Complex> other, Matrix<Complex> result)
        {
            // diagonal - diagonal = diagonal
            var diagOther = other as DiagonalMatrix;
            var diagResult = result as DiagonalMatrix;
            if (diagOther != null && diagResult != null)
            {
                Control.LinearAlgebraProvider.SubtractArrays(_data, diagOther._data, diagResult._data);
                return;
            }

            other.Negate(result);
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, result.At(i, i) + _data[i - 1]);
            }
        }

        /// <summary>
        /// Multiplies each element of the matrix by a scalar and places results into the result matrix.
        /// </summary>
        /// <param name="scalar">The scalar to multiply the matrix with.</param>
        /// <param name="result">The matrix to store the result of the multiplication.</param>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        protected override void DoMultiply(Complex scalar, Matrix<Complex> result)
        {
            if (scalar == 0.0)
            {
                result.Clear();
                return;
            }

            if (scalar == 1.0)
            {
                CopyTo(result);
                return;
            }

            var diagResult = result as DiagonalMatrix;
            if (diagResult == null)
            {
                base.DoMultiply(scalar, result);
            }
            else
            {
                Control.LinearAlgebraProvider.ScaleArray(scalar, _data, diagResult._data);
            }
        }

        /// <summary>
        /// Multiplies this matrix with a vector and places the results into the result vector.
        /// </summary>
        /// <param name="rightSide">The vector to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoMultiply(Vector<Complex> rightSide, Vector<Complex> result)
        {
            var d = Math.Min(ColumnCount, RowCount);
            if (d < RowCount)
            {
                result.ClearSubVector(ColumnCount, RowCount - ColumnCount);
            }

            if (d == ColumnCount)
            {
                var denseOther = rightSide.Storage as DenseVectorStorage<Complex>;
                var denseResult = result.Storage as DenseVectorStorage<Complex>;
                if (denseOther != null && denseResult != null)
                {
                    Control.LinearAlgebraProvider.PointWiseMultiplyArrays(_data, denseOther.Data, denseResult.Data);
                    return;
                }
            }

            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, _data[i - 1]*rightSide.At(i));
            }
        }

        /// <summary>
        /// Multiplies this matrix with another matrix and places the results into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoMultiply(Matrix<Complex> other, Matrix<Complex> result)
        {
            var diagonalOther = other as DiagonalMatrix;
            var diagonalResult = result as DiagonalMatrix;
            if (diagonalOther != null && diagonalResult != null)
            {
                var thisDataCopy = new Complex[diagonalResult._data.Length];
                var otherDataCopy = new Complex[diagonalResult._data.Length];
                Array.Copy(_data, thisDataCopy, (diagonalResult._data.Length > _data.Length) ? _data.Length : diagonalResult._data.Length);
                Array.Copy(diagonalOther._data, otherDataCopy, (diagonalResult._data.Length > diagonalOther._data.Length) ? diagonalOther._data.Length : diagonalResult._data.Length);
                Control.LinearAlgebraProvider.PointWiseMultiplyArrays(thisDataCopy, otherDataCopy, diagonalResult._data);
                return;
            }

            var denseOther = other.Storage as DenseColumnMajorMatrixStorage<Complex>;
            if (denseOther != null)
            {
                var dense = denseOther.Data;
                var diagonal = _data;
                var d = Math.Min(denseOther.RowCount, RowCount);
                if (d < RowCount)
                {
                    result.ClearSubMatrix(denseOther.RowCount, RowCount - denseOther.RowCount, 1, denseOther.ColumnCount);
                }
                int index = 0;
                for (int i = 1; i <= denseOther.ColumnCount; i++)
                {
                    for (int j = 1; j <= d; j++)
                    {
                        result.At(j, i, dense[index]*diagonal[j - 1]);
                        index++;
                    }
                    index += (denseOther.RowCount - d);
                }
                return;
            }

            if (ColumnCount == RowCount)
            {
                other.Storage.MapIndexedTo(result.Storage, (i, j, x) => x*_data[i - 1], Zeros.AllowSkip, ExistingData.Clear);
            }
            else
            {
                result.Clear();
                other.Storage.MapSubMatrixIndexedTo(result.Storage, (i, j, x) => x*_data[i - 1], 1, 1, other.RowCount, 1, 1, other.ColumnCount, Zeros.AllowSkip, ExistingData.AssumeZeros);
            }
        }

        /// <summary>
        /// Multiplies this matrix with transpose of another matrix and places the results into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoTransposeAndMultiply(Matrix<Complex> other, Matrix<Complex> result)
        {
            var diagonalOther = other as DiagonalMatrix;
            var diagonalResult = result as DiagonalMatrix;
            if (diagonalOther != null && diagonalResult != null)
            {
                var thisDataCopy = new Complex[diagonalResult._data.Length];
                var otherDataCopy = new Complex[diagonalResult._data.Length];
                Array.Copy(_data, thisDataCopy, (diagonalResult._data.Length > _data.Length) ? _data.Length : diagonalResult._data.Length);
                Array.Copy(diagonalOther._data, otherDataCopy, (diagonalResult._data.Length > diagonalOther._data.Length) ? diagonalOther._data.Length : diagonalResult._data.Length);
                Control.LinearAlgebraProvider.PointWiseMultiplyArrays(thisDataCopy, otherDataCopy, diagonalResult._data);
                return;
            }

            var denseOther = other.Storage as DenseColumnMajorMatrixStorage<Complex>;
            if (denseOther != null)
            {
                var dense = denseOther.Data;
                var diagonal = _data;
                var d = Math.Min(denseOther.ColumnCount, RowCount);
                if (d < RowCount)
                {
                    result.ClearSubMatrix(denseOther.ColumnCount, RowCount - denseOther.ColumnCount, 1, denseOther.RowCount);
                }
                int index = 0;
                for (int j = 1; j <= d; j++)
                {
                    var diagj = diagonal[j - 1];
                    for (int i = 1; i <= denseOther.RowCount; i++)
                    {
                        result.At(j, i, dense[index]*diagj);
                        index++;
                    }
                }
                return;
            }

            base.DoTransposeAndMultiply(other, result);
        }

        /// <summary>
        /// Multiplies this matrix with the conjugate transpose of another matrix and places the results into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoConjugateTransposeAndMultiply(Matrix<Complex> other, Matrix<Complex> result)
        {
            var diagonalOther = other as DiagonalMatrix;
            var diagonalResult = result as DiagonalMatrix;
            if (diagonalOther != null && diagonalResult != null)
            {
                var thisDataCopy = new Complex[diagonalResult._data.Length];
                var otherDataCopy = new Complex[diagonalResult._data.Length];
                Array.Copy(_data, thisDataCopy, (diagonalResult._data.Length > _data.Length) ? _data.Length : diagonalResult._data.Length);
                Array.Copy(diagonalOther._data, otherDataCopy, (diagonalResult._data.Length > diagonalOther._data.Length) ? diagonalOther._data.Length : diagonalResult._data.Length);
                // TODO: merge/MulByConj
                Control.LinearAlgebraProvider.ConjugateArray(otherDataCopy, otherDataCopy);
                Control.LinearAlgebraProvider.PointWiseMultiplyArrays(thisDataCopy, otherDataCopy, diagonalResult._data);
                return;
            }

            var denseOther = other.Storage as DenseColumnMajorMatrixStorage<Complex>;
            if (denseOther != null)
            {
                var dense = denseOther.Data;
                var diagonal = _data;
                var d = Math.Min(denseOther.ColumnCount, RowCount);
                if (d < RowCount)
                {
                    result.ClearSubMatrix(denseOther.ColumnCount, RowCount - denseOther.ColumnCount, 1, denseOther.RowCount);
                }
                int index = 0;
                for (int j = 1; j <= d; j++)
                {
                    var diagj = diagonal[j - 1];
                    for (int i = 1; i <= denseOther.RowCount; i++)
                    {
                        result.At(j, i, dense[index].Conjugate()*diagj);
                        index++;
                    }
                }
                return;
            }

            base.DoConjugateTransposeAndMultiply(other, result);
        }

        /// <summary>
        /// Multiplies the transpose of this matrix with another matrix and places the results into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoTransposeThisAndMultiply(Matrix<Complex> other, Matrix<Complex> result)
        {
            var diagonalOther = other as DiagonalMatrix;
            var diagonalResult = result as DiagonalMatrix;
            if (diagonalOther != null && diagonalResult != null)
            {
                var thisDataCopy = new Complex[diagonalResult._data.Length];
                var otherDataCopy = new Complex[diagonalResult._data.Length];
                Array.Copy(_data, thisDataCopy, (diagonalResult._data.Length > _data.Length) ? _data.Length : diagonalResult._data.Length);
                Array.Copy(diagonalOther._data, otherDataCopy, (diagonalResult._data.Length > diagonalOther._data.Length) ? diagonalOther._data.Length : diagonalResult._data.Length);
                Control.LinearAlgebraProvider.PointWiseMultiplyArrays(thisDataCopy, otherDataCopy, diagonalResult._data);
                return;
            }

            var denseOther = other.Storage as DenseColumnMajorMatrixStorage<Complex>;
            if (denseOther != null)
            {
                var dense = denseOther.Data;
                var diagonal = _data;
                var d = Math.Min(denseOther.RowCount, ColumnCount);
                if (d < ColumnCount)
                {
                    result.ClearSubMatrix(denseOther.RowCount, ColumnCount - denseOther.RowCount, 1, denseOther.ColumnCount);
                }
                int index = 0;
                for (int i = 1; i <= denseOther.ColumnCount; i++)
                {
                    for (int j = 1; j <= d; j++)
                    {
                        result.At(j, i, dense[index] * diagonal[j - 1]);
                        index++;
                    }
                    index += (denseOther.RowCount - d);
                }
                return;
            }

            if (ColumnCount == RowCount)
            {
                other.Storage.MapIndexedTo(result.Storage, (i, j, x) => x*_data[i - 1], Zeros.AllowSkip, ExistingData.Clear);
            }
            else
            {
                result.Clear();
                other.Storage.MapSubMatrixIndexedTo(result.Storage, (i, j, x) => x*_data[i - 1], 1, 1, other.RowCount, 1, 1, other.ColumnCount, Zeros.AllowSkip, ExistingData.AssumeZeros);
            }
        }

        /// <summary>
        /// Multiplies the transpose of this matrix with another matrix and places the results into the result matrix.
        /// </summary>
        /// <param name="other">The matrix to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoConjugateTransposeThisAndMultiply(Matrix<Complex> other, Matrix<Complex> result)
        {
            var diagonalOther = other as DiagonalMatrix;
            var diagonalResult = result as DiagonalMatrix;
            if (diagonalOther != null && diagonalResult != null)
            {
                var thisDataCopy = new Complex[diagonalResult._data.Length];
                var otherDataCopy = new Complex[diagonalResult._data.Length];
                Array.Copy(_data, thisDataCopy, (diagonalResult._data.Length > _data.Length) ? _data.Length : diagonalResult._data.Length);
                Array.Copy(diagonalOther._data, otherDataCopy, (diagonalResult._data.Length > diagonalOther._data.Length) ? diagonalOther._data.Length : diagonalResult._data.Length);
                // TODO: merge/MulByConj
                Control.LinearAlgebraProvider.ConjugateArray(thisDataCopy, thisDataCopy);
                Control.LinearAlgebraProvider.PointWiseMultiplyArrays(thisDataCopy, otherDataCopy, diagonalResult._data);
                return;
            }

            var denseOther = other.Storage as DenseColumnMajorMatrixStorage<Complex>;
            if (denseOther != null)
            {
                var dense = denseOther.Data;
                var conjugateDiagonal = new Complex[_data.Length];
                for (int i = 0; i < _data.Length; i++)
                {
                    conjugateDiagonal[i] = _data[i].Conjugate();
                }

                var d = Math.Min(denseOther.RowCount, ColumnCount);
                if (d < ColumnCount)
                {
                    result.ClearSubMatrix(denseOther.RowCount, ColumnCount - denseOther.RowCount, 1, denseOther.ColumnCount);
                }
                int index = 0;
                for (int i = 1; i <= denseOther.ColumnCount; i++)
                {
                    for (int j = 1; j <= d; j++)
                    {
                        result.At(j, i, dense[index]*conjugateDiagonal[j - 1]);
                        index++;
                    }
                    index += (denseOther.RowCount - d);
                }
                return;
            }

            base.DoConjugateTransposeThisAndMultiply(other, result);
        }

        /// <summary>
        /// Multiplies the transpose of this matrix with a vector and places the results into the result vector.
        /// </summary>
        /// <param name="rightSide">The vector to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoTransposeThisAndMultiply(Vector<Complex> rightSide, Vector<Complex> result)
        {
            var d = Math.Min(ColumnCount, RowCount);
            if (d < ColumnCount)
            {
                result.ClearSubVector(RowCount, ColumnCount - RowCount);
            }

            if (d == RowCount)
            {
                var denseOther = rightSide.Storage as DenseVectorStorage<Complex>;
                var denseResult = result.Storage as DenseVectorStorage<Complex>;
                if (denseOther != null && denseResult != null)
                {
                    Control.LinearAlgebraProvider.PointWiseMultiplyArrays(_data, denseOther.Data, denseResult.Data);
                    return;
                }
            }

            for (var i = 1; i <= d; i++)
            {
                result.At(i, _data[i - 1]*rightSide.At(i));
            }
        }

        /// <summary>
        /// Multiplies the conjugate transpose of this matrix with a vector and places the results into the result vector.
        /// </summary>
        /// <param name="rightSide">The vector to multiply with.</param>
        /// <param name="result">The result of the multiplication.</param>
        protected override void DoConjugateTransposeThisAndMultiply(Vector<Complex> rightSide, Vector<Complex> result)
        {
            var d = Math.Min(ColumnCount, RowCount);
            if (d < ColumnCount)
            {
                result.ClearSubVector(RowCount, ColumnCount - RowCount);
            }

            if (d == RowCount)
            {
                var denseOther = rightSide.Storage as DenseVectorStorage<Complex>;
                var denseResult = result.Storage as DenseVectorStorage<Complex>;
                if (denseOther != null && denseResult != null)
                {
                    // TODO: merge/MulByConj
                    Control.LinearAlgebraProvider.ConjugateArray(_data, denseResult.Data);
                    Control.LinearAlgebraProvider.PointWiseMultiplyArrays(denseResult.Data, denseOther.Data, denseResult.Data);
                    return;
                }
            }

            for (var i = 1; i <= d; i++)
            {
                result.At(i, _data[i - 1].Conjugate()*rightSide.At(i));
            }
        }

        /// <summary>
        /// Divides each element of the matrix by a scalar and places results into the result matrix.
        /// </summary>
        /// <param name="divisor">The scalar to divide the matrix with.</param>
        /// <param name="result">The matrix to store the result of the division.</param>
        protected override void DoDivide(Complex divisor, Matrix<Complex> result)
        {
            if (divisor == Complex.One)
            {
                CopyTo(result);
                return;
            }

            var diagResult = result as DiagonalMatrix;
            if (diagResult != null)
            {
                Control.LinearAlgebraProvider.ScaleArray(1.0/divisor, _data, diagResult._data);
                return;
            }

            result.Clear();
            for (int i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, _data[i - 1]/divisor);
            }
        }

        /// <summary>
        /// Divides a scalar by each element of the matrix and stores the result in the result matrix.
        /// </summary>
        /// <param name="dividend">The scalar to add.</param>
        /// <param name="result">The matrix to store the result of the division.</param>
        protected override void DoDivideByThis(Complex dividend, Matrix<Complex> result)
        {
            var diagResult = result as DiagonalMatrix;
            if (diagResult != null)
            {
                var resultData = diagResult._data;
                CommonParallel.For(0, _data.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        resultData[i] = dividend/_data[i];
                    }
                });
                return;
            }

            result.Clear();
            for (int i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, dividend/_data[i - 1]);
            }
        }

        /// <summary>
        /// Computes the determinant of this matrix.
        /// </summary>
        /// <returns>The determinant of this matrix.</returns>
        public override Complex Determinant()
        {
            if (RowCount != ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            return _data.Aggregate(Complex.One, (current, t) => current * t);
        }

        /// <summary>
        /// Returns the elements of the diagonal in a <see cref="DenseVector"/>.
        /// </summary>
        /// <returns>The elements of the diagonal.</returns>
        /// <remarks>For non-square matrices, the method returns Min(Rows, Columns) elements where
        /// i == j (i is the row index, and j is the column index).</remarks>
        public override Vector<Complex> Diagonal()
        {
            return new DenseVector(_data).Clone();
        }

        /// <summary>
        /// Copies the values of the given array to the diagonal.
        /// </summary>
        /// <param name="source">The array to copy the values from. The length of the vector should be
        /// Min(Rows, Columns).</param>
        /// <exception cref="ArgumentException">If the length of <paramref name="source"/> does not
        /// equal Min(Rows, Columns).</exception>
        /// <remarks>For non-square matrices, the elements of <paramref name="source"/> are copied to
        /// this[i,i].</remarks>
        public override void SetDiagonal(Complex[] source)
        {
            if (source.Length != _data.Length)
            {
                throw new ArgumentException(Resources.ArgumentArraysSameLength, "source");
            }

            Array.Copy(source, _data, source.Length);
        }

        /// <summary>
        /// Copies the values of the given <see cref="Vector{T}"/> to the diagonal.
        /// </summary>
        /// <param name="source">The vector to copy the values from. The length of the vector should be
        /// Min(Rows, Columns).</param>
        /// <exception cref="ArgumentException">If the length of <paramref name="source"/> does not
        /// equal Min(Rows, Columns).</exception>
        /// <remarks>For non-square matrices, the elements of <paramref name="source"/> are copied to
        /// this[i,i].</remarks>
        public override void SetDiagonal(Vector<Complex> source)
        {
            var denseSource = source as DenseVector;
            if (denseSource == null)
            {
                base.SetDiagonal(source);
                return;
            }

            if (_data.Length != denseSource.Values.Length)
            {
                throw new ArgumentException(Resources.ArgumentVectorsSameLength, "source");
            }

            Array.Copy(denseSource.Values, _data, denseSource.Values.Length);
        }

        /// <summary>Calculates the induced L1 norm of this matrix.</summary>
        /// <returns>The maximum absolute column sum of the matrix.</returns>
        public override double L1Norm()
        {
            return _data.Aggregate(0d, (current, t) => Math.Max(current, t.Magnitude));
        }

        /// <summary>Calculates the induced L2 norm of the matrix.</summary>
        /// <returns>The largest singular value of the matrix.</returns>
        public override double L2Norm()
        {
            return _data.Aggregate(0d, (current, t) => Math.Max(current, t.Magnitude));
        }

        /// <summary>Calculates the induced infinity norm of this matrix.</summary>
        /// <returns>The maximum absolute row sum of the matrix.</returns>
        public override double InfinityNorm()
        {
            return L1Norm();
        }

        /// <summary>Calculates the Frobenius norm of this matrix.</summary>
        /// <returns>The Frobenius norm of this matrix.</returns>
        public override double FrobeniusNorm()
        {
            return Math.Sqrt(_data.Sum(t => t.Magnitude * t.Magnitude));
        }

        /// <summary>Calculates the condition number of this matrix.</summary>
        /// <returns>The condition number of the matrix.</returns>
        public override Complex ConditionNumber()
        {
            var maxSv = double.NegativeInfinity;
            var minSv = double.PositiveInfinity;
            foreach (var t in _data)
            {
                maxSv = Math.Max(maxSv, t.Magnitude);
                minSv = Math.Min(minSv, t.Magnitude);
            }

            return maxSv / minSv;
        }

        /// <summary>Computes the inverse of this matrix.</summary>
        /// <exception cref="ArgumentException">If <see cref="DiagonalMatrix"/> is not a square matrix.</exception>
        /// <exception cref="ArgumentException">If <see cref="DiagonalMatrix"/> is singular.</exception>
        /// <returns>The inverse of this matrix.</returns>
        public override Matrix<Complex> Inverse()
        {
            if (RowCount != ColumnCount)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSquare);
            }

            var inverse = (DiagonalMatrix)Clone();
            for (var i = 0; i < _data.Length; i++)
            {
                if (_data[i] != 0.0)
                {
                    inverse._data[i] = 1.0 / _data[i];
                }
                else
                {
                    throw new ArgumentException(Resources.ArgumentMatrixNotSingular);
                }
            }

            return inverse;
        }

        /// <summary>
        /// Returns a new matrix containing the lower triangle of this matrix.
        /// </summary>
        /// <returns>The lower triangle of this matrix.</returns>
        public override Matrix<Complex> LowerTriangle()
        {
            return Clone();
        }

        /// <summary>
        /// Puts the lower triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public override void LowerTriangle(Matrix<Complex> result)
        {
            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw DimensionsDontMatch<ArgumentException>(this, result, "result");
            }

            if (ReferenceEquals(this, result))
            {
                return;
            }

            result.Clear();
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, _data[i - 1]);
            }
        }

        /// <summary>
        /// Returns a new matrix containing the lower triangle of this matrix. The new matrix
        /// does not contain the diagonal elements of this matrix.
        /// </summary>
        /// <returns>The lower triangle of this matrix.</returns>
        public override Matrix<Complex> StrictlyLowerTriangle()
        {
            return new DiagonalMatrix(RowCount, ColumnCount);
        }

        /// <summary>
        /// Puts the strictly lower triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public override void StrictlyLowerTriangle(Matrix<Complex> result)
        {
            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw DimensionsDontMatch<ArgumentException>(this, result, "result");
            }

            result.Clear();
        }

        /// <summary>
        /// Returns a new matrix containing the upper triangle of this matrix.
        /// </summary>
        /// <returns>The upper triangle of this matrix.</returns>
        public override Matrix<Complex> UpperTriangle()
        {
            return Clone();
        }

        /// <summary>
        /// Puts the upper triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public override void UpperTriangle(Matrix<Complex> result)
        {
            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw DimensionsDontMatch<ArgumentException>(this, result, "result");
            }

            result.Clear();
            for (var i = 1; i <= _data.Length; i++)
            {
                result.At(i, i, _data[i - 1]);
            }
        }

        /// <summary>
        /// Returns a new matrix containing the upper triangle of this matrix. The new matrix
        /// does not contain the diagonal elements of this matrix.
        /// </summary>
        /// <returns>The upper triangle of this matrix.</returns>
        public override Matrix<Complex> StrictlyUpperTriangle()
        {
            return new DiagonalMatrix(RowCount, ColumnCount);
        }

        /// <summary>
        /// Puts the strictly upper triangle of this matrix into the result matrix.
        /// </summary>
        /// <param name="result">Where to store the lower triangle.</param>
        /// <exception cref="ArgumentException">If the result matrix's dimensions are not the same as this matrix.</exception>
        public override void StrictlyUpperTriangle(Matrix<Complex> result)
        {
            if (result.RowCount != RowCount || result.ColumnCount != ColumnCount)
            {
                throw DimensionsDontMatch<ArgumentException>(this, result, "result");
            }

            result.Clear();
        }

        /// <summary>
        /// Creates a matrix that contains the values from the requested sub-matrix.
        /// </summary>
        /// <param name="rowIndex">The row to start copying from.</param>
        /// <param name="rowCount">The number of rows to copy. Must be positive.</param>
        /// <param name="columnIndex">The column to start copying from.</param>
        /// <param name="columnCount">The number of columns to copy. Must be positive.</param>
        /// <returns>The requested sub-matrix.</returns>
        /// <exception cref="ArgumentOutOfRangeException">If: <list><item><paramref name="rowIndex"/> is
        /// &lt; 1, or greater than the number of rows.</item>
        /// <item><paramref name="columnIndex"/> is &lt; 1, or greater than the number
        /// of columns.</item>
        /// <item><c>(columnIndex + columnLength - 1) &gt; Columns</c></item>
        /// <item><c>(rowIndex + rowLength - 1) &gt; Rows</c></item></list></exception>
        /// <exception cref="ArgumentOutOfRangeException">If <paramref name="rowCount"/> or <paramref name="columnCount"/>
        /// is negative.</exception>
        public override Matrix<Complex> SubMatrix(int rowIndex, int rowCount, int columnIndex, int columnCount)
        {
            var target = rowIndex == columnIndex
                ? (Matrix<Complex>)new DiagonalMatrix(rowCount, columnCount)
                : new SparseMatrix(rowCount, columnCount);

            Storage.CopySubMatrixTo(target.Storage, rowIndex, 1, rowCount, columnIndex, 1, columnCount, ExistingData.AssumeZeros);
            return target;
        }

        /// <summary>
        /// Permute the columns of a matrix according to a permutation.
        /// </summary>
        /// <param name="p">The column permutation to apply to this matrix.</param>
        /// <exception cref="InvalidOperationException">Always thrown</exception>
        /// <remarks>Permutation in diagonal matrix are senseless, because of matrix nature</remarks>
        public override void PermuteColumns(Permutation p)
        {
            throw new InvalidOperationException("Permutations in diagonal matrix are not allowed");
        }

        /// <summary>
        /// Permute the rows of a matrix according to a permutation.
        /// </summary>
        /// <param name="p">The row permutation to apply to this matrix.</param>
        /// <exception cref="InvalidOperationException">Always thrown</exception>
        /// <remarks>Permutation in diagonal matrix are senseless, because of matrix nature</remarks>
        public override void PermuteRows(Permutation p)
        {
            throw new InvalidOperationException("Permutations in diagonal matrix are not allowed");
        }

        /// <summary>
        /// Evaluates whether this matrix is symmetric.
        /// </summary>
        public override sealed bool IsSymmetric()
        {
            return true;
        }

        /// <summary>
        /// Evaluates whether this matrix is hermitian (conjugate symmetric).
        /// </summary>
        public override sealed bool IsHermitian()
        {
            for (var k = 0; k < _data.Length; k ++)
            {
                if (!_data[k].IsReal())
                {
                    return false;
                }
            }

            return true;
        }
    }
}