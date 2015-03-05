﻿// <copyright file="DenseVectorStorage.cs" company="Math.NET">
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
using System.Linq;
using MathNet.Numerics.Properties;
using MathNet.Numerics.Threading;

namespace MathNet.Numerics.LinearAlgebra.OneBased.Storage
{
    [Serializable]
    public class DenseVectorStorage<T> : VectorStorage<T>
        where T : struct, IEquatable<T>, IFormattable
    {
        // [ruegg] public fields are OK here

        public readonly T[] Data;

        internal DenseVectorStorage(int length)
            : base(length)
        {
            Data = new T[length];
        }

        internal DenseVectorStorage(int length, T[] data)
            : base(length)
        {
            if (data == null)
            {
                throw new ArgumentNullException("data");
            }

            if (data.Length != length)
            {
                throw new ArgumentOutOfRangeException("data", string.Format(Resources.ArgumentArrayWrongLength, length));
            }

            Data = data;
        }

        /// <summary>
        /// True if the vector storage format is dense.
        /// </summary>
        public override bool IsDense
        {
            get { return true; }
        }

        /// <summary>
        /// Retrieves the requested element without range checking.
        /// </summary>
        public override T At(int index)
        {
            return Data[index - 1];
        }

        /// <summary>
        /// Sets the element without range checking.
        /// </summary>
        public override void At(int index, T value)
        {
            Data[index - 1] = value;
        }

        public override void Clear()
        {
            Array.Clear(Data, 0, Data.Length);
        }

        public override void Clear(int index, int count)
        {
            Array.Clear(Data, index - 1, count);
        }

        // INITIALIZATION

        public static DenseVectorStorage<T> OfVector(VectorStorage<T> vector)
        {
            var storage = new DenseVectorStorage<T>(vector.Length);
            vector.CopyToUnchecked(storage, ExistingData.AssumeZeros);
            return storage;
        }

        public static DenseVectorStorage<T> OfValue(int length, T value)
        {
            if (length < 0)
            {
                throw new ArgumentOutOfRangeException("length", Resources.ArgumentNotNegative);
            }

            var data = new T[length];
            CommonParallel.For(0, data.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    data[i] = value;
                }
            });
            return new DenseVectorStorage<T>(length, data);
        }

        public static DenseVectorStorage<T> OfInit(int length, Func<int, T> init)
        {
            if (length < 0)
            {
                throw new ArgumentOutOfRangeException("length", Resources.ArgumentNotNegative);
            }

            if (init == null)
            {
                throw new ArgumentNullException("init");
            }

            var data = new T[length];
            CommonParallel.For(0, data.Length, 4096, (a, b) =>
            {
                for (int i = a; i < b; i++)
                {
                    data[i] = init(i + 1);
                }
            });
            return new DenseVectorStorage<T>(length, data);
        }

        public static DenseVectorStorage<T> OfEnumerable(IEnumerable<T> data)
        {
            if (data == null)
            {
                throw new ArgumentNullException("data");
            }

            var arrayData = data as T[];
            if (arrayData != null)
            {
                var copy = new T[arrayData.Length];
                Array.Copy(arrayData, copy, arrayData.Length);
                return new DenseVectorStorage<T>(copy.Length, copy);
            }

            var array = data.ToArray();
            return new DenseVectorStorage<T>(array.Length, array);
        }

        public static DenseVectorStorage<T> OfIndexedEnumerable(int length, IEnumerable<Tuple<int, T>> data)
        {
            if (data == null)
            {
                throw new ArgumentNullException("data");
            }

            var array = new T[length];
            foreach (var item in data)
            {
                array[item.Item1 - 1] = item.Item2;
            }
            return new DenseVectorStorage<T>(array.Length, array);
        }

        // VECTOR COPY

        internal override void CopyToUnchecked(VectorStorage<T> target, ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseVectorStorage<T>;
            if (denseTarget != null)
            {
                if (!ReferenceEquals(this, denseTarget))
                {
                    Array.Copy(Data, 0, denseTarget.Data, 0, Data.Length);
                }

                return;
            }

            var sparseTarget = target as SparseVectorStorage<T>;
            if (sparseTarget != null)
            {
                var indices = new List<int>();
                var values = new List<T>();

                for (int i = 0; i < Data.Length; i++)
                {
                    var item = Data[i];
                    if (!Zero.Equals(item))
                    {
                        values.Add(item);
                        indices.Add(i);
                    }
                }

                sparseTarget.Indices = indices.ToArray();
                sparseTarget.Values = values.ToArray();
                sparseTarget.ValueCount = values.Count;
                return;
            }

            // FALL BACK

            for (int i = 0; i < Data.Length; i++)
            {
                target.At(i + 1, Data[i]);
            }
        }

        // ROW COPY

        internal override void CopyToRowUnchecked(MatrixStorage<T> target, int rowIndex, ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseColumnMajorMatrixStorage<T>;
            if (denseTarget != null)
            {
                --rowIndex;     // adjust for use in 0-based indexing calculation
                for (int j = 0; j < Data.Length; j++)
                {
                    denseTarget.Data[j*target.RowCount + rowIndex] = Data[j];
                }
                return;
            }

            // FALL BACK

            for (int j = 0; j < Length; j++)
            {
                target.At(rowIndex, j, Data[j]);
            }
        }

        // COLUMN COPY

        internal override void CopyToColumnUnchecked(MatrixStorage<T> target, int columnIndex, ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseColumnMajorMatrixStorage<T>;
            if (denseTarget != null)
            {
                // adjust columnIndex for use in 0-based indexing calculation
                Array.Copy(Data, 0, denseTarget.Data, (columnIndex - 1)*denseTarget.RowCount, Data.Length);
                return;
            }

            // FALL BACK

            for (int i = 0; i < Length; i++)
            {
                target.At(i + 1, columnIndex, Data[i]);
            }
        }

        // SUB-VECTOR COPY

        internal override void CopySubVectorToUnchecked(VectorStorage<T> target,
            int sourceIndex, int targetIndex, int count,
            ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseVectorStorage<T>;
            if (denseTarget != null)
            {
                // adjust sourceIndex and targetIndex for use in 0-based indexing calculation
                Array.Copy(Data, sourceIndex - 1, denseTarget.Data, targetIndex - 1, count);
                return;
            }

            // FALL BACK

            base.CopySubVectorToUnchecked(target, sourceIndex, targetIndex, count, existingData);
        }

        // SUB-ROW COPY

        internal override void CopyToSubRowUnchecked(MatrixStorage<T> target, int rowIndex,
            int sourceColumnIndex, int targetColumnIndex, int columnCount,
            ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseColumnMajorMatrixStorage<T>;
            if (denseTarget != null)
            {
                // adjust sourceColumnIndex, targetColumnIndex and rowIndex for use in 0-based indexing calculation
                --sourceColumnIndex;
                --targetColumnIndex;
                --rowIndex;
                for (int j = 0; j < Data.Length; j++)
                {
                    denseTarget.Data[(j + targetColumnIndex)*target.RowCount + rowIndex] = Data[j + sourceColumnIndex];
                }
                return;
            }

            // FALL BACK

            // adjust sourceColumnIndex for use in 0-based indexing calculation
            --sourceColumnIndex;
            for (int j = sourceColumnIndex, jj = targetColumnIndex; j < sourceColumnIndex + columnCount; j++, jj++)
            {
                target.At(rowIndex, jj, Data[j]);
            }
        }

        // SUB-COLUMN COPY

        internal override void CopyToSubColumnUnchecked(MatrixStorage<T> target, int columnIndex,
            int sourceRowIndex, int targetRowIndex, int rowCount,
            ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseColumnMajorMatrixStorage<T>;
            if (denseTarget != null)
            {
                Array.Copy(Data, sourceRowIndex - 1, denseTarget.Data, (columnIndex - 1)*denseTarget.RowCount + targetRowIndex - 1, rowCount);
                return;
            }

            // FALL BACK

            // adjust sourceRowIndex for use in 0-based indexing calculation
            --sourceRowIndex;
            for (int i = sourceRowIndex, ii = targetRowIndex; i < sourceRowIndex + rowCount; i++, ii++)
            {
                target.At(ii, columnIndex, Data[i]);
            }
        }

        // ENUMERATION

        public override IEnumerable<T> Enumerate()
        {
            return Data;
        }

        public override IEnumerable<Tuple<int, T>> EnumerateIndexed()
        {
            return Data.Select((t, i) => new Tuple<int, T>(i + 1, t));
        }

        public override IEnumerable<T> EnumerateNonZero()
        {
            return Data.Where(x => !Zero.Equals(x));
        }

        public override IEnumerable<Tuple<int, T>> EnumerateNonZeroIndexed()
        {
            for (var i = 0; i < Data.Length; i++)
            {
                if (!Zero.Equals(Data[i]))
                {
                    yield return new Tuple<int, T>(i + 1, Data[i]);
                }
            }
        }

        // FUNCTIONAL COMBINATORS

        internal override void MapToUnchecked<TU>(VectorStorage<TU> target, Func<T, TU> f,
            Zeros zeros = Zeros.AllowSkip, ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseVectorStorage<TU>;
            if (denseTarget != null)
            {
                CommonParallel.For(0, Data.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        denseTarget.Data[i] = f(Data[i]);
                    }
                });
                return;
            }

            // FALL BACK

            for (int i = 0; i < Length; i++)
            {
                target.At(i + 1, f(Data[i]));
            }
        }

        internal override void MapIndexedToUnchecked<TU>(VectorStorage<TU> target, Func<int, T, TU> f,
            Zeros zeros = Zeros.AllowSkip, ExistingData existingData = ExistingData.Clear)
        {
            var denseTarget = target as DenseVectorStorage<TU>;
            if (denseTarget != null)
            {
                CommonParallel.For(0, Data.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        denseTarget.Data[i] = f(i + 1, Data[i]);
                    }
                });
                return;
            }

            // FALL BACK

            for (int i = 0, i1 = 1; i < Length; i++, i1++)
            {
                target.At(i1, f(i1, Data[i]));
            }
        }

        internal override void Map2ToUnchecked(VectorStorage<T> target, VectorStorage<T> other, Func<T, T, T> f, Zeros zeros = Zeros.AllowSkip, ExistingData existingData = ExistingData.Clear)
        {
            if (target is SparseVectorStorage<T>)
            {
                // Recursive to dense target at first, since the operation is
                // effectively dense anyway because at least one operand is dense
                var intermediate = new DenseVectorStorage<T>(target.Length);
                Map2ToUnchecked(intermediate, other, f, zeros, ExistingData.AssumeZeros);
                intermediate.CopyTo(target, existingData);
                return;
            }

            var denseTarget = target as DenseVectorStorage<T>;
            var denseOther = other as DenseVectorStorage<T>;
            if (denseTarget != null && denseOther != null)
            {
                CommonParallel.For(0, Data.Length, 4096, (a, b) =>
                {
                    for (int i = a; i < b; i++)
                    {
                        denseTarget.Data[i] = f(Data[i], denseOther.Data[i]);
                    }
                });

                return;
            }

            var sparseOther = other as SparseVectorStorage<T>;
            if (denseTarget != null && sparseOther != null)
            {
                T[] targetData = denseTarget.Data;
                int[] otherIndices = sparseOther.Indices;
                T[] otherValues = sparseOther.Values;
                int otherValueCount = sparseOther.ValueCount;

                int k = 0;
                for (int i = 0; i < Data.Length; i++)
                {
                    if (k < otherValueCount && otherIndices[k] == i)
                    {
                        targetData[i] = f(Data[i], otherValues[k]);
                        k++;
                    }
                    else
                    {
                        targetData[i] = f(Data[i], Zero);
                    }
                }

                return;
            }

            base.Map2ToUnchecked(target, other, f, zeros, existingData);
        }

        internal override TState Fold2Unchecked<TOther, TState>(VectorStorage<TOther> other, Func<TState, T, TOther, TState> f, TState state, Zeros zeros = Zeros.AllowSkip)
        {
            var denseOther = other as DenseVectorStorage<TOther>;
            if (denseOther != null)
            {
                var otherData = denseOther.Data;
                for (int i = 0; i < Data.Length; i++)
                {
                    state = f(state, Data[i], otherData[i]);
                }

                return state;
            }

            var sparseOther = other as SparseVectorStorage<TOther>;
            if (sparseOther != null)
            {
                int[] otherIndices = sparseOther.Indices;
                TOther[] otherValues = sparseOther.Values;
                int otherValueCount = sparseOther.ValueCount;
                TOther otherZero = BuilderInstance<TOther>.Vector.Zero;

                int k = 0;
                for (int i = 0; i < Data.Length; i++)
                {
                    if (k < otherValueCount && otherIndices[k] == i)
                    {
                        state = f(state, Data[i], otherValues[k]);
                        k++;
                    }
                    else
                    {
                        state = f(state, Data[i], otherZero);
                    }
                }

                return state;
            }

            return base.Fold2Unchecked(other, f, state, zeros);
        }
    }
}
