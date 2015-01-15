// <copyright file="VectorStorage.Validation.cs" company="Math.NET">
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

namespace MathNet.Numerics.LinearAlgebra.OneBased.Storage
{
    // ReSharper disable UnusedParameter.Global
    public partial class VectorStorage<T>
    {
        void ValidateRange(int index)
        {
            if (index < 1 || index > Length)
            {
                throw new ArgumentOutOfRangeException("index");
            }
        }

        void ValidateSubVectorRange(VectorStorage<T> target,
            int sourceIndex, int targetIndex, int count)
        {
            if (count < 0)
            {
                throw new ArgumentOutOfRangeException("count", Resources.ArgumentNotNegative);
            }

            // Verify Source

            if (sourceIndex > Length || sourceIndex < 1)
            {
                throw new ArgumentOutOfRangeException("sourceIndex");
            }

            var sourceMax = sourceIndex + count - 1;

            if (sourceMax > Length)
            {
                throw new ArgumentOutOfRangeException("count");
            }

            // Verify Target

            if (targetIndex > target.Length || targetIndex < 1)
            {
                throw new ArgumentOutOfRangeException("targetIndex");
            }

            var targetMax = targetIndex + count - 1;

            if (targetMax > target.Length)
            {
                throw new ArgumentOutOfRangeException("count");
            }
        }

        void ValidateRowRange(MatrixStorage<T> target, int rowIndex)
        {
            if (rowIndex > target.RowCount || rowIndex < 1)
            {
                throw new ArgumentOutOfRangeException("rowIndex");
            }

            if (target.ColumnCount != Length)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameRowDimension, "target");
            }
        }

        void ValidateColumnRange(MatrixStorage<T> target, int columnIndex)
        {
            if (columnIndex > target.ColumnCount || columnIndex < 1)
            {
                throw new ArgumentOutOfRangeException("columnIndex");
            }

            if (target.RowCount != Length)
            {
                throw new ArgumentException(Resources.ArgumentMatrixSameColumnDimension, "target");
            }
        }

        void ValidateSubRowRange(MatrixStorage<T> target, int rowIndex,
            int sourceColumnIndex, int targetColumnIndex, int columnCount)
        {
            if (columnCount < 0)
            {
                throw new ArgumentOutOfRangeException("columnCount", Resources.ArgumentMustBePositive);
            }

            // Verify Source

            if (sourceColumnIndex > Length || sourceColumnIndex < 1)
            {
                throw new ArgumentOutOfRangeException("sourceColumnIndex");
            }

            if (sourceColumnIndex + columnCount - 1 > Length)
            {
                throw new ArgumentOutOfRangeException("columnCount");
            }

            // Verify Target

            if (rowIndex > target.RowCount - 1 || rowIndex < 1)
            {
                throw new ArgumentOutOfRangeException("rowIndex");
            }

            if (targetColumnIndex > target.ColumnCount || targetColumnIndex < 1)
            {
                throw new ArgumentOutOfRangeException("targetColumnIndex");
            }

            if (targetColumnIndex + columnCount - 1 > target.ColumnCount)
            {
                throw new ArgumentOutOfRangeException("columnCount");
            }
        }

        void ValidateSubColumnRange(MatrixStorage<T> target, int columnIndex,
            int sourceRowIndex, int targetRowIndex, int rowCount)
        {
            if (rowCount < 0)
            {
                throw new ArgumentOutOfRangeException("rowCount", Resources.ArgumentMustBePositive);
            }

            // Verify Source

            if (sourceRowIndex > Length || sourceRowIndex < 1)
            {
                throw new ArgumentOutOfRangeException("sourceRowIndex");
            }

            if (sourceRowIndex + rowCount - 1 > Length)
            {
                throw new ArgumentOutOfRangeException("rowCount");
            }

            // Verify Target

            if (columnIndex > target.ColumnCount || columnIndex < 1)
            {
                throw new ArgumentOutOfRangeException("columnIndex");
            }

            if (targetRowIndex > target.RowCount || targetRowIndex < 1)
            {
                throw new ArgumentOutOfRangeException("targetRowIndex");
            }

            if (targetRowIndex + rowCount - 1 > target.RowCount)
            {
                throw new ArgumentOutOfRangeException("rowCount");
            }
        }
    }
    // ReSharper restore UnusedParameter.Global
}
