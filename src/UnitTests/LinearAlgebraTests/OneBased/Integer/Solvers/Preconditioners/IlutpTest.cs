// <copyright file="IlutpTest.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// http://mathnetnumerics.codeplex.com
//
// Copyright (c) 2009-2013 Math.NET
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
using MathNet.Numerics.LinearAlgebra.Integer.Solvers;
using NUnit.Framework;

namespace MathNet.Numerics.UnitTests.LinearAlgebraTests.OneBased.Integer.Solvers.Preconditioners
{
    /// <summary>
    /// Incomplete LU with IlutpPreconditioner test with drop tolerance and partial pivoting.
    /// </summary>
    [TestFixture, Category("LASolver")]
    public sealed class IlutpPreconditionerTest
    {
        /// <summary>
        /// The drop tolerance.
        /// </summary>
        double _dropTolerance = 0.1;

        /// <summary>
        /// The fill level.
        /// </summary>
        double _fillLevel = 1.0;

        /// <summary>
        /// The pivot tolerance.
        /// </summary>
        double _pivotTolerance = 1.0;

        /// <summary>
        /// Default constructor throws <c>NotSupportedException</c>.
        /// </summary>
        [Test]
        public void DefaultConstructorThrowsNotSupportedException()
        {
            Assert.Throws<NotSupportedException>(() => {
                var preconditioner = new ILUTPPreconditioner();
            });
        }

        /// <summary>
        /// Parameterized constructor throws <c>NotSupportedException</c>.
        /// </summary>
        [Test]
        public void ParameterizedConstructorThrowsNotSupportedException()
        {
            Assert.Throws<NotSupportedException>(() => {
                var preconditioner = new ILUTPPreconditioner(_fillLevel,_dropTolerance,_pivotTolerance);
            });
        }
    }
}
