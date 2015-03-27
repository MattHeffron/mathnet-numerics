using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Solvers;
using MathNet.Numerics.LinearAlgebra.Complex.Solvers;
using MathNet.Numerics.LinearAlgebra.Complex;
using System.Numerics;
using System.Diagnostics;

namespace ConsoleApplication0
{
    class Program
    {
        static void Main(string[] args)
        {
            ////for (int i = 0; i < 15; i++)
            ////{
            ////    var vec = Vector<int>.Build.DenseOfEnumerable(Enumerable.Range(1, i));
            ////    string foo = vec.ToVectorString(3,5);
            ////}
            SetupMatrices();
            //CanMultiplyMatrixWithMatrix("Tall3x2", "Wide2x3");
            CanComputeDeterminant("Square3x3");
        }

        /// <summary>
        /// Gets or sets test matrices values to use.
        /// </summary>
        static Dictionary<string, Complex[,]> TestData2D { get; set; }

        /// <summary>
        /// Gets or sets test matrices instances to use.
        /// </summary>
        static Dictionary<string, Matrix<Complex>> TestMatrices { get; set; }

        public static void SetupMatrices()
        {
            TestData2D = new Dictionary<string, Complex[,]>
                {
                    {"Singular3x3", new[,] {{new Complex(1.0, 1), Complex.Zero, Complex.Zero}, {Complex.Zero, Complex.Zero, Complex.Zero}, {Complex.Zero, Complex.Zero, new Complex(3.0, 1)}}},
                    {"Square3x3", new[,] {{new Complex(-1.1, 1), Complex.Zero, Complex.Zero}, {Complex.Zero, new Complex(1.1, 1), Complex.Zero}, {Complex.Zero, Complex.Zero, new Complex(6.6, 1)}}},
                    {"Square4x4", new[,] {{new Complex(-1.1, 1), Complex.Zero, Complex.Zero, Complex.Zero}, {Complex.Zero, new Complex(1.1, 1), Complex.Zero, Complex.Zero}, {Complex.Zero, Complex.Zero, new Complex(6.2, 1), Complex.Zero}, {Complex.Zero, Complex.Zero, Complex.Zero, new Complex(-7.7, 1)}}},
                    {"Singular4x4", new[,] {{new Complex(-1.1, 1), Complex.Zero, Complex.Zero, Complex.Zero}, {Complex.Zero, new Complex(-2.2, 1), Complex.Zero, Complex.Zero}, {Complex.Zero, Complex.Zero, Complex.Zero, Complex.Zero}, {Complex.Zero, Complex.Zero, Complex.Zero, new Complex(-4.4, 1)}}},
                    {"Tall3x2", new[,] {{new Complex(-1.1, 1), Complex.Zero}, {Complex.Zero, new Complex(1.1, 1)}, {Complex.Zero, Complex.Zero}}},
                    {"Wide2x3", new[,] {{new Complex(-1.1, 1), Complex.Zero, Complex.Zero}, {Complex.Zero, new Complex(1.1, 1), Complex.Zero}}}
                };

            TestMatrices = new Dictionary<string, Matrix<Complex>>();
            foreach (var name in TestData2D.Keys)
            {
                TestMatrices.Add(name, DiagonalMatrix.OfArray(TestData2D[name]));
            }
        }

        public static void CanComputeDeterminant(string name)
        {
            var matrix = TestMatrices[name];
            var denseMatrix = DenseMatrix.OfArray(TestData2D[name]);
            Debug.WriteLine("matrix: " + matrix.ToString());
            Debug.WriteLine("denseMatrix: " + denseMatrix.ToString());
            var dmDet = denseMatrix.Determinant();
            var mDet = matrix.Determinant();
        }
        public static void CanMultiplyMatrixWithMatrix(string nameA, string nameB)
        {
            var matrixA = TestMatrices[nameA];
            Debug.WriteLine("matrixA: " + matrixA.ToString());
            var matrixB = TestMatrices[nameB];
            Debug.WriteLine("matrixB: " + matrixB.ToString());
            var matrixC = matrixA * matrixB;
            Debug.WriteLine("matrixC: " + matrixC.ToString());

            Debug.Assert(matrixC.RowCount == matrixA.RowCount);
            Debug.Assert(matrixC.ColumnCount == matrixB.ColumnCount);
            for (int i = 1; i <= matrixC.RowCount; i++)
            {
                for (int j = 1; j <= matrixC.ColumnCount; j++)
                {
                    var raw = matrixA.Row(i) * matrixB.Column(j);
                    var dot = matrixA.Row(i).DotProduct(matrixB.Column(j));
                    var mC = matrixC[i, j];
                    var diff = mC - raw;
                }
            }
        }
    }
}
