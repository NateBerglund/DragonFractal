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
        /// <param name="color">Base color to use for rendering</param>
        public void RenderFractal(DirectBitmap image, int nIterations, Action<Matrix<double>, DirectBitmap, int> renderUnit, Matrix<double> transform, int color)
        {
            if (nIterations > 0)
                foreach (Matrix<double> fracFunc in IFS)
                    RenderFractal(image, nIterations - 1, renderUnit, transform * fracFunc, color);
            else
                renderUnit(transform, image, color);
        }

        /// <summary>
        /// Draws a line from (-1, 1) to (1, 1) (transformed according to transform) as the base unit.
        /// </summary>
        /// <param name="transform">Coordinate transform to apply</param>
        /// <param name="image">Image upon which to render the line</param>
        /// <param name="color">Color to use</param>
        public static void renderLine(Matrix<double> transform, DirectBitmap image, int color)
        {
            Vector<double> pointA = DenseVector.OfArray(new double[] { -1.0, 1.0, 1.0 });
            Vector<double> pointB = DenseVector.OfArray(new double[] { 1.0, 1.0, 1.0 });
            Vector<double> pointAT = transform * pointA;
            Vector<double> pointBT = transform * pointB;
            ImProc.DrawLine(pointAT[0] / pointAT[2], pointAT[1] / pointAT[2], pointBT[0] / pointBT[2], pointBT[1] / pointBT[2], color, image);
        }

        /// <summary>
        /// Renders the primary spirals of the fractal in red (transformed according to transform) as the base unit.
        /// NOTE: This will only render the spirals correctly if the transform is a product of
        /// rotations, translations and scalings (must scale both coordinates uniformly)
        /// </summary>
        /// <param name="transform">Coordinate transform to apply</param>
        /// <param name="image">Image upon which to render the spirals</param>
        /// <param name="color">Color to use</param>
        public static void renderPrimarySpirals(Matrix<double> transform, DirectBitmap image, int color)
        {
            Vector<double> pointA = DenseVector.OfArray(new double[] { -1.0, 1.0, 1.0 });
            Vector<double> pointB = DenseVector.OfArray(new double[] { 1.0, 1.0, 1.0 });
            Vector<double> pointAT = transform * pointA;
            Vector<double> pointBT = transform * pointB;
            double scaleFactor = Math.Sqrt(transform.SubMatrix(0, 2, 0, 2).Determinant());
            Vector<double> xBasis = DenseVector.OfArray(new double[] { 1.0, 0.0, 0.0 }); // make 3rd coordinate zero so translation is ignored
            Vector<double> xBasisT = transform * xBasis;
            double theta = Math.Atan2(xBasisT[1], xBasisT[0]);
            ImProc.DrawSpiral(pointAT[0] / pointAT[2], pointAT[1] / pointAT[2], theta + Math.PI / 2, scaleFactor * 1.0 / 3.0, 6 * Math.PI, Math.PI / 2, Math.PI / 512, 1.0 / 16.0, color, image);
            ImProc.DrawSpiral(pointBT[0] / pointBT[2], pointBT[1] / pointBT[2], theta - Math.PI / 2, scaleFactor * 2.0 / 3.0, 6 * Math.PI, Math.PI / 2, Math.PI / 512, 1.0 / 16.0, color, image);
        }

        /// <summary>
        /// Renders the secondary spirals of the fractal in blue (transformed according to transform) as the base unit.
        /// NOTE: This will only render the spirals correctly if the transform is a product of
        /// rotations, translations and scalings (must scale both coordinates uniformly)
        /// </summary>
        /// <param name="transform">Coordinate transform to apply</param>
        /// <param name="image">Image upon which to render the spirals</param>
        /// <param name="color">Color to use</param>
        /// <param name="color">Color to use</param>
        public static void renderSecondarySpirals(Matrix<double> transform, DirectBitmap image, int color)
        {
            Vector<double> origin = DenseVector.OfArray(new double[] { 0.0, 0.0, 1.0 });
            Vector<double> originT = transform * origin;
            double scaleFactor = Math.Sqrt(transform.SubMatrix(0, 2, 0, 2).Determinant());
            Vector<double> xBasis = DenseVector.OfArray(new double[] { 1.0, 0.0, 0.0 }); // make 3rd coordinate zero so translation is ignored
            Vector<double> xBasisT = transform * xBasis;
            double theta = Math.Atan2(xBasisT[1], xBasisT[0]);
            ImProc.DrawSpiral(originT[0] / originT[2], originT[1] / originT[2], theta + Math.PI / 2, scaleFactor * 1.0 / 2.0, 6 * Math.PI, Math.PI / 4, Math.PI / 512, 1.0 / 16.0, color, image);
        }
    }
}
