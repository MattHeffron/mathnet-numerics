using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra.OneBased;
using MathNet.Numerics.LinearAlgebra.OneBased.Solvers;
using MathNet.Numerics.LinearAlgebra.OneBased.Complex.Solvers;
using MathNet.Numerics.LinearAlgebra.OneBased.Complex;
using System.Numerics;
using System.Diagnostics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;

namespace ConsoleApplication1
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
            //CanNormalizeRows(1);
            //CanMultiplyMatrixWithMatrix("Tall3x2", "Wide2x3");
            //CanComputeDeterminant("Square3x3");
            CanFactorizeRandomMatrix(1);
        }
        public static void CanFactorizeRandomMatrix(int order)
        {
            var matrixX = Matrix<Complex>.Build.RandomPositiveDefinite(order, 1);
            var chol = matrixX.Cholesky();
            var factorC = chol.Factor;

            // Make sure the Cholesky factor has the right dimensions.
            //Assert.AreEqual(order, factorC.RowCount);
            //Assert.AreEqual(order, factorC.ColumnCount);

            //// Make sure the Cholesky factor is lower triangular.
            //AssertHelpers.IsLowerTriangular(factorC);

            //// Make sure the cholesky factor times it's transpose is the original matrix.
            //var matrixXfromC = factorC * factorC.ConjugateTranspose();
            //AssertHelpers.AlmostEqualRelative(matrixX, matrixXfromC, 8);
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
        public static void CanNormalizeRows(int p)
        {
            var matrix = TestMatrices["Square4x4"].NormalizeRows(p);
            //for (var i = 1; i <= matrix.RowCount; i++)
            //{
            //    var row = matrix.Row(i);
            //    AssertHelpers.AlmostEqual(Complex.One, row.Norm(p), 12);
            //}
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
