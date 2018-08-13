﻿using System;
using System.Collections.Generic;
using System.IO;
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
            string OUTPUT_DIRECTORY = @"C:\Users\Dianna\data\DragonFractal";
            Directory.CreateDirectory(OUTPUT_DIRECTORY);
            double base_size = (512.0 / 30.0); // base size is approximately 17.067 inches, fractal will be approximately 3 * base_size by 2 * base_size.
            double subImageWidth = 9.0; // width of sub-images in inches
            double subImageHeight = 6.5; // height of sub-images in inches
            int previewSF = 16; // Scale factor to use between the initial low-res preview and the final rull-res result. Must be a power of 2.

            // Iterated function system defining the dragon fractal
            Fractal fractal = new Fractal();
            fractal.IFS = new List<Matrix<double>>();
            var Translate = DenseMatrix.OfArray(new double[,] { { 1, 0, 1 }, { 0, 1, -1 }, { 0, 0, 1 } });
            var Scale = DenseMatrix.OfArray(new double[,] { { Math.Sqrt(0.5), 0, 0 }, { 0, Math.Sqrt(0.5), 0 }, { 0, 0, 1 } });
            var Rot45 = DenseMatrix.OfArray(new double[,] { { Math.Cos(Math.PI / 4), -Math.Sin(Math.PI / 4), 0 }, { Math.Sin(Math.PI / 4), Math.Cos(Math.PI / 4), 0 }, { 0, 0, 1 } });
            var Rot135 = DenseMatrix.OfArray(new double[,] { { Math.Cos(3 * Math.PI / 4), -Math.Sin(3 * Math.PI / 4), 0 }, { Math.Sin(3 * Math.PI / 4), Math.Cos(3 * Math.PI / 4), 0 }, { 0, 0, 1 } });
            fractal.IFS.Add(Rot45 * Scale * Translate);
            fractal.IFS.Add(Rot135 * Scale * Translate);

            double[] xCoords = null, yCoords = null;

            for (int bigIter = 0; bigIter < 2; ++bigIter)
            {
                int sf = (0 == bigIter ? 1 : previewSF) * 256; // note: needs to be a power of 2
                int DPI = (int)Math.Round(sf / base_size); // dots (pixels) per inch
                int pad = (0 == bigIter ? 1 : previewSF) * 140;
                int w = 3 * sf + pad;
                int h = 2 * sf + pad;
                var finalTransform = DenseMatrix.OfArray(new double[,] { { sf, 0, 4 * sf / 3 + pad / 2 }, { 0, sf, sf / 3 + pad / 2 }, { 0, 0, 1 } });
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

                    // Render fractal to 18 iterations, using a simple line segment as the base unit
                    fractal.RenderFractal(image, 0 == bigIter ? 18 : 26, Fractal.renderLine, finalTransform, black, black);
                    ImProc.BWBoundary(image, borderImage, black); // compute boundary
                    fractal.ReferenceImage = borderImage;
                    fractal.AuxImage = thickImage;

                    // Render the primary spirals of the fractal onto both image and thickImage
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

                    int nIterSS = 0 == bigIter ? 10 : 14; // number of iterations of the secondary spirals to compute
                    // Render iterations of the secondary spirals onto the fractal
                    for (int i = 0; i < nIterSS; ++i)
                        fractal.RenderFractal(image, i, fractal.renderSecondarySpirals, finalTransform, blue, black);

                    // Dilate the first image with a disc shaped kernel of radius 1/32 inch
                    // Then draw these pixels in black onto the thickImage. This will make the spirals of the thick image
                    // never get thinner than our smallest tool.
                    {
                        byte[,] tempMask = ImProc.ToBWByteEq(image, blue);
                        bool[,] kernel = ImProc.GenerateDiscKernel(DPI / 32.0);
                        int kernelRadius = (kernel.GetLength(0) - 1) / 2;
                        tempMask = ImProc.BWDilate(tempMask, kernel, kernelRadius, kernelRadius);
                        ImProc.DrawColor(thickImage, tempMask, black);
                    }
                    GC.Collect(); // free up memory we temporarily allocated

                    // Clean up the white diagonal lines
                    ImProc.BWFilter(thickImage, ImProc.CleanWhiteDiagonalsLUT);

                    IO.SaveImage(image, Path.Combine(OUTPUT_DIRECTORY, "thinfractal" + bigIter.ToString() + ".bmp"));
                    IO.SaveImage(thickImage, Path.Combine(OUTPUT_DIRECTORY, "fractal" + bigIter.ToString() + ".bmp"));
                    IO.SaveImage(borderImage, Path.Combine(OUTPUT_DIRECTORY, "fractalborder" + bigIter.ToString() + ".bmp"));

                    if (0 == bigIter)
                    {
                        // Dilate fractal with a circular filter of radius 60
                        byte[,] maskImage = ImProc.ToBWByteNeq(image, white); // mask indicates where original image is non-white
                        bool[,] kernel = ImProc.GenerateDiscKernel(60);
                        int kernelRadius = (kernel.GetLength(0) - 1) / 2;
                        maskImage = ImProc.BWDilate(maskImage, kernel, kernelRadius, kernelRadius);
                        byte[,] maskImage2 = new byte[h, w];
                        ImProc.BWBoundary4(maskImage, maskImage2, 255);
                        maskImage = null;
                        // Detect the contour
                        ImProc.FindClosedContour(maskImage2, out xCoords, out yCoords);
                        // Smooth the contour
                        double sigma = 5.0;
                        int halfWidth = (int)Math.Ceiling(3.0 * sigma);
                        double[] gaussKernel = ImProc.ConstructGaussianKernel(sigma, halfWidth, true);
                        double[] xCoordsSmooth = ImProc.Convolve1D(xCoords, gaussKernel, 1, -1, 1);
                        double[] yCoordsSmooth = ImProc.Convolve1D(yCoords, gaussKernel, 1, -1, 1);
                        // Decimate smoothed curves so curvature based smoothing will work faster
                        int N = xCoords.Length;
                        int N2 = N / 10;
                        xCoords = new double[N2];
                        yCoords = new double[N2];
                        for (int i = 0; i < N2; ++i)
                        {
                            xCoords[i] = xCoordsSmooth[(int)(i * N / (double)N2 + 0.5)];
                            yCoords[i] = yCoordsSmooth[(int)(i * N / (double)N2 + 0.5)];
                        }
                        // Apply further curvature-based smoothing
                        ImProc.CurvatureSmoothContour(xCoords, yCoords, 2000.0, 1.0, 10.0, true);
                    }
                    else
                    {
                        // Scale contour by a factor of previewSF
                        int N = xCoords.Length;
                        for (int i = 0; i < N; ++i)
                        {
                            xCoords[i] = xCoords[i] * previewSF;
                            yCoords[i] = yCoords[i] * previewSF;
                        }
                    }

                    // Extract boundary of thick fractal
                    ImProc.SetColor(image, black);
                    ImProc.BWBoundary(thickImage, image, black);

                    // Open the fractal with a disc shaped kernel of radius 1/16 inch,
                    // adding this boundary to the original boundary. This will show
                    // us how far a tool of radius 1/16 inch can make it into the spirals.
                    {
                        byte[,] tempMask = ImProc.ToBWByteEq(thickImage, black);
                        bool[,] kernel = ImProc.GenerateDiscKernel(DPI / 16.0);
                        int kernelRadius = (kernel.GetLength(0) - 1) / 2;
                        tempMask = ImProc.BWOpen(tempMask, kernel, kernelRadius, kernelRadius);
                        ImProc.BWBoundary(tempMask, image, 255);
                    }
                    GC.Collect(); // free up memory we temporarily allocated

                    // Open the fractal with a disc shaped kernel of radius 1/8 inch,
                    // adding this boundary to the original boundary. This will show
                    // us how far a tool of radius 1/8 inch can make it into the spirals.
                    {
                        byte[,] tempMask = ImProc.ToBWByteEq(thickImage, black);
                        bool[,] kernel = ImProc.GenerateDiscKernel(DPI / 8.0);
                        int kernelRadius = (kernel.GetLength(0) - 1) / 2;
                        tempMask = ImProc.BWOpen(tempMask, kernel, kernelRadius, kernelRadius);
                        ImProc.BWBoundary(tempMask, image, 255);
                    }
                    GC.Collect(); // free up memory we temporarily allocated

                    // Open the fractal with a disc shaped kernel of radius 1/4 inch,
                    // adding this boundary to the original boundary. This will show
                    // us how far a tool of radius 1/4 inch can make it into the spirals.
                    {
                        byte[,] tempMask = ImProc.ToBWByteEq(thickImage, black);
                        bool[,] kernel = ImProc.GenerateDiscKernel(DPI / 4.0);
                        int kernelRadius = (kernel.GetLength(0) - 1) / 2;
                        tempMask = ImProc.BWOpen(tempMask, kernel, kernelRadius, kernelRadius);
                        ImProc.BWBoundary(tempMask, image, 255);
                    }
                    GC.Collect(); // free up memory we temporarily allocated

                    // Invert the image
                    ImProc.BWInvert(image);

                    // Draw contour onto image
                    ImProc.DrawClosedContour(xCoords, yCoords, black, image);

                    //// Add fractal boundary to image
                    //for (int y = 0; y < h; ++y)
                    //    for (int x = 0; x < w; ++x)
                    //        if (borderImage.Bits[w * y + x] == white)
                    //            image.Bits[w * y + x] = black;

                    IO.SaveImage(image, Path.Combine(OUTPUT_DIRECTORY, "fractalwithboundary" + bigIter.ToString() + ".bmp"));

                    // Divide the image into sub-images and save them
                    if (1 == bigIter)
                    {
                        DirectBitmap[] subImages = ImProc.partitionImage(image, (int)Math.Round(subImageWidth * DPI), (int)Math.Round(subImageHeight * DPI), black);
                        int nImages = subImages.Length;
                        for (int i = 0; i < nImages; ++i)
                            IO.SaveImage(subImages[i], Path.Combine(OUTPUT_DIRECTORY, "subimage" + i.ToString() + ".bmp"), DPI, DPI);
                    }
                }
            }
        }
    }
}
