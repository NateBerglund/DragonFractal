using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace DragonFractal
{
    class Program
    {
        static void Main(string[] args)
        {
            Fractal fractal = new Fractal();
            fractal.IFS = new List<Matrix<double>>();
            var Translate = DenseMatrix.OfArray(new double[,] { { 1, 0, 1 }, { 0, 1, -1 }, { 0, 0, 1 } });
            var Scale = DenseMatrix.OfArray(new double[,] { { Math.Sqrt(0.5), 0, 0 }, { 0, Math.Sqrt(0.5), 0 }, { 0, 0, 1 } });
            var Rot45 = DenseMatrix.OfArray(new double[,] { { Math.Cos(Math.PI/4), -Math.Sin(Math.PI / 4), 0 }, { Math.Sin(Math.PI / 4), Math.Cos(Math.PI / 4), 0 }, { 0, 0, 1 } });
            var Rot135 = DenseMatrix.OfArray(new double[,] { { Math.Cos(3 * Math.PI / 4), -Math.Sin(3 * Math.PI / 4), 0 }, { Math.Sin(3 * Math.PI / 4), Math.Cos(3 * Math.PI / 4), 0 }, { 0, 0, 1 } });
            fractal.IFS.Add(Rot45 * Scale * Translate);
            fractal.IFS.Add(Rot135 * Scale * Translate);
            int sf = 256; // 2048
            int w = 3 * sf + 100;
            int h = 2 * sf + 100;
            var finalTransform = DenseMatrix.OfArray(new double[,] { { sf, 0, 4 * sf / 3 + 50 }, { 0, sf, sf / 3 + 50 }, { 0, 0, 1 } });
            using (DirectBitmap image = new DirectBitmap(w, h),
                borderImage = new DirectBitmap(w, h),
                thickImage = new DirectBitmap(w, h))
            {
                int white = unchecked((int)0xffffffff); // white
                int black = unchecked((int)0xff000000); // black
                int red = unchecked((int)0xffff0000); // red
                int blue = unchecked((int)0xff0000ff); // blue
                //int green = unchecked((int)0xff00ff00); // green

                // Initialize the image with white
                for (int y = 0; y < h; ++y)
                {
                    for (int x = 0; x < w; ++x)
                    {
                        image.Bits[w * y + x] = white;
                        thickImage.Bits[w * y + x] = white;
                    }
                }

                fractal.RenderFractal(image, 18, Fractal.renderLine, finalTransform, black, black);
                ImProc.BWBoundary(image, borderImage, black);
                fractal.ReferenceImage = borderImage;
                fractal.AuxImage = thickImage;

                //fractal.RenderFractal(image, 17, Fractal.renderLine, finalTransform * fractal.IFS[0], green);
                fractal.RenderFractal(image, 0, Fractal.renderPrimarySpirals, finalTransform, blue, black);
                fractal.RenderFractal(thickImage, 0, Fractal.renderPrimarySpirals, finalTransform, black, black);

                // Draw special thick spirals directly onto the thickImage
                Vector<double> pointA = DenseVector.OfArray(new double[] { -1.0, 1.0, 1.0 });
                Vector<double> pointB = DenseVector.OfArray(new double[] { 1.0, 1.0, 1.0 });
                Vector<double> pointAT = finalTransform * pointA;
                Vector<double> pointBT = finalTransform * pointB;
                double scaleFactor = Math.Sqrt(finalTransform.SubMatrix(0, 2, 0, 2).Determinant());
                Vector<double> xBasis = DenseVector.OfArray(new double[] { 1.0, 0.0, 0.0 }); // make 3rd coordinate zero so translation is ignored
                Vector<double> xBasisT = finalTransform * xBasis;
                double theta = Math.Atan2(xBasisT[1], xBasisT[0]);
                for (int i = 0; i < 12; ++i)
                {
                    ImProc.DrawTaperedSpiral(pointAT[0] / pointAT[2], pointAT[1] / pointAT[2], theta + Math.PI / 4 + i * Math.PI / 4, scaleFactor * Math.Sqrt(2.0) * Math.Pow(0.5, 0.5 * i) / 3.0, 0, Math.PI / 4, Math.PI / 1000, 1.0 / 16.0, 0.02, 1.0, black, thickImage);
                    ImProc.DrawTaperedSpiral(pointBT[0] / pointBT[2], pointBT[1] / pointBT[2], theta - 3 * Math.PI / 4 + i * Math.PI / 4, scaleFactor * Math.Sqrt(2.0) * Math.Pow(0.5, 0.5 * i) * 2.0 / 3.0, 0, Math.PI / 4, Math.PI / 1000, 1.0 / 16.0, 0.02, 1.0, black, thickImage);
                }

                for (int i = 0; i < 10; ++i)
                    fractal.RenderFractal(image, i, fractal.renderSecondarySpirals, finalTransform, blue, black);

                IO.SaveImage(image, @"C:\Users\info\Desktop\thinfractal.bmp");
                IO.SaveImage(thickImage, @"C:\Users\info\Desktop\fractal.bmp");
                IO.SaveImage(borderImage, @"C:\Users\info\Desktop\fractalborder.bmp");

                // Dilate fractal with a circular filter of radius 40
                byte[,] maskImage = new byte[h, w];
                for (int y = 0; y < h; ++y)
                    for (int x = 0; x < w; ++x)
                        if (image.Bits[w * y + x] != white)
                            maskImage[y, x] = 255;
                int kernelRadius = 40;
                bool[,] kernel = new bool[2 * kernelRadius + 1, 2 * kernelRadius + 1];
                for (int y = -kernelRadius; y <= kernelRadius; ++y)
                    for (int x = -kernelRadius; x <= kernelRadius; ++x)
                        if (x * x + y * y <= kernelRadius * kernelRadius)
                            kernel[kernelRadius + y, kernelRadius + x] = true;
                maskImage = ImProc.BWDilate(maskImage, kernel, kernelRadius, kernelRadius);
                byte[,] maskImage2 = new byte[h, w];
                ImProc.BWBoundary4(maskImage, maskImage2, 255);
                maskImage = null;
                // Detect the contour
                double[] xCoords, yCoords;
                ImProc.FindClosedContour(maskImage2, out xCoords, out yCoords);
                // Smooth the contour
                double sigma = 5.0;
                int halfWidth = (int)Math.Ceiling(3.0 * sigma);
                double[] gaussKernel = ImProc.ConstructGaussianKernel(sigma, halfWidth, true);
                xCoords = ImProc.Convolve1D(xCoords, gaussKernel, 1, -1, 1);
                yCoords = ImProc.Convolve1D(yCoords, gaussKernel, 1, -1, 1);
                // Apply further curvature-based smoothing


                // Draw contour onto thickImage
                for (int y = 0; y < h; ++y)
                    for (int x = 0; x < w; ++x)
                        image.Bits[w * y + x] = black;
                ImProc.DrawClosedContour(xCoords, yCoords, black, thickImage);

                IO.SaveImage(thickImage, @"C:\Users\info\Desktop\fractalwithboundary.bmp");
            }
        }
    }
}
