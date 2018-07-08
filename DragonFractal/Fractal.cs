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
        /// Image to use for reference. Can be used by base unit drawing functions
        /// to selectively suppress drawing.
        /// </summary>
        public DirectBitmap ReferenceImage = null;

        /// <summary>
        /// Auxillary image to use for drawing. Can be used by base unit drawing functions
        /// to selectively suppress drawing.
        /// </summary>
        public DirectBitmap AuxImage = null;

        /// <summary>
        /// Fractal rendering function
        /// </summary>
        /// <param name="image">Image to render fractal onto</param>
        /// <param name="nIterations">Number of iterations to perform</param>
        /// <param name="renderUnit">Function that renders a single basic unit that IFS will be applied to recursively</param>
        /// <param name="transform">Matrix that transforms from the fractal's coordinate system to the coordinate system of the image</param>
        /// <param name="color">Base color to use for rendering</param>
        /// <param name="secondaryColor">Color to render onto auxilliary image (not used)</param>
        public void RenderFractal(DirectBitmap image, int nIterations, Action<Matrix<double>, DirectBitmap, int, int> renderUnit, Matrix<double> transform, int color, int secondaryColor)
        {
            if (nIterations > 0)
            {
                for (int i = 0; i < IFS.Count; ++i)
                {
                    Matrix<double> fracFunc = IFS[i];
                    RenderFractal(image, nIterations - 1, renderUnit, transform * fracFunc, color, secondaryColor);
                }
            }
            else
                renderUnit(transform, image, color, secondaryColor);
        }

        /// <summary>
        /// Draws a line from (-1, 1) to (1, 1) (transformed according to transform) as the base unit.
        /// </summary>
        /// <param name="transform">Coordinate transform to apply</param>
        /// <param name="image">Image upon which to render the line</param>
        /// <param name="color">Color to use</param>
        /// <param name="secondaryColor">Color to render onto auxilliary image (not used)</param>
        public static void renderLine(Matrix<double> transform, DirectBitmap image, int color, int secondaryColor)
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
        /// <param name="secondaryColor">Color to render onto auxilliary image (not used)</param>
        public static void renderPrimarySpirals(Matrix<double> transform, DirectBitmap image, int color, int secondaryColor)
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
        /// <param name="secondaryColor">Color to render onto auxilliary image</param>
        public void renderSecondarySpirals(Matrix<double> transform, DirectBitmap image, int color, int secondaryColor)
        {    
            Vector<double> origin = DenseVector.OfArray(new double[] { 0.0, 0.0, 1.0 });
            Vector<double> originT = transform * origin;
            // NOTE: Second spiral originating at pointB is not strictly necessary
            // (it will be constructed as a copy in another iteration). We only include
            // it so we can extend the theta bounds farther than they would otherwise go).
            Vector<double> pointB = DenseVector.OfArray(new double[] { -0.5, 0.5, 1.0 });
            Vector<double> pointBT = transform * pointB;
            double scaleFactor = Math.Sqrt(transform.SubMatrix(0, 2, 0, 2).Determinant());
            Vector<double> xBasis = DenseVector.OfArray(new double[] { 1.0, 0.0, 0.0 }); // make 3rd coordinate zero so translation is ignored
            Vector<double> xBasisT = transform * xBasis;
            double theta = Math.Atan2(xBasisT[1], xBasisT[0]);
            // Suppress rendering if the origin is not a boundary pixel of the dragon fractal
            bool render = true;
            if (null != ReferenceImage)
            {
                int x = (int)(originT[0] / originT[2] + 0.5);
                int y = (int)(originT[1] / originT[2] + 0.5);
                if (x >= 0 && x < ReferenceImage.Width && y >= 0 && y < ReferenceImage.Height)
                {
                    if ((ReferenceImage.Bits[ReferenceImage.Width * y + x] & 0x00ffffff) == 0)
                        render = false; 
                }
            }
            if (render)
            {
                double thetaPlusSpan, thetaMinusSpan;
                ImProc.DrawSpiral(originT[0] / originT[2], originT[1] / originT[2], theta + Math.PI, scaleFactor * 1.0 / 4.0, Math.PI / 512, 1.0 / 16.0, color, image, out thetaPlusSpan, out thetaMinusSpan);
                if (null != AuxImage && (thetaMinusSpan > 0 || thetaPlusSpan > 0))
                    ImProc.DrawThickSpiral(originT[0] / originT[2], originT[1] / originT[2], theta + Math.PI, scaleFactor * 1.0 / 4.0, thetaPlusSpan, thetaMinusSpan, Math.PI / 512, 1.0 / 16.0, 0.1, secondaryColor, AuxImage);
            }
        }
    }
}
