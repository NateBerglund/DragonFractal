using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace DragonFractal
{
    /// <summary>
    /// Class used for generating fractals
    /// </summary>
    public class Fractal
    {
        /// <summary>
        /// Iterated function system defining the fractal. This is a list of 3x3
        /// matrices in homogenous coordinates.
        /// </summary>
        public List<Matrix<double>> IFS;

        /// <summary>
        /// Fractal rendering function
        /// </summary>
        /// <param name="image">Image to render fractal onto</param>
        /// <param name="nIterations">Number of iterations to perform</param>
        /// <param name="renderUnit">Function that renders a single basic unit that IFS will be applied to recursively</param>
        /// <param name="transform">Matrix that transforms from the fractal's coordinate system to the coordinate system of the image</param>
        public void RenderFractal(DirectBitmap image, int nIterations, Action<Matrix<double>, DirectBitmap> renderUnit, Matrix<double> transform)
        {
            if (nIterations > 0)
                foreach (Matrix<double> fracFunc in IFS)
                    RenderFractal(image, nIterations - 1, renderUnit, transform * fracFunc);
            else
                renderUnit(transform, image);
        }

        /// <summary>
        /// Draws a black line from (-1, 1) to (1, 1) (transformed according to transform) as the base unit.
        /// </summary>
        /// <param name="transform">Coordinate transform to apply</param>
        /// <param name="image">Image upon which to render the line</param>
        public static void renderLine(Matrix<double> transform, DirectBitmap image)
        {
            int color = unchecked((int)0xff000000); // black
            Vector<double> pointA = DenseVector.OfArray(new double[] { -1.0, 1.0, 1.0 });
            Vector<double> pointB = DenseVector.OfArray(new double[] { 1.0, 1.0, 1.0 });
            Vector<double> pointAT = transform * pointA;
            Vector<double> pointBT = transform * pointB;
            ImProc.DrawLine(pointAT[0] / pointAT[2], pointAT[1] / pointAT[2], pointBT[0] / pointBT[2], pointBT[1] / pointBT[2], color, image);
        }
    }
}
