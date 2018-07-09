using System;
using System.Collections.Generic;

namespace DragonFractal
{
    /// <summary>
    /// Class for performing various image processing operations
    /// </summary>
    static class ImProc
    {
        #region Drawing

        /// <summary>
        /// Draws a single point on an image. No subpixel sampling is used.
        /// </summary>
        /// <param name="x">x-coordinate of the point to draw</param>
        /// <param name="y">y-coordinate of the point to draw</param>
        /// <param name="color">Color to give the point</param>
        /// <param name="image">Image to draw the point onto</param>
        public static void DrawPt(double x, double y, int color, DirectBitmap image)
        {
            int xi = (int)(x + 0.5);
            int yi = (int)(y + 0.5);
            if (xi >= 0 && xi < image.Width && yi >= 0 && yi < image.Height)
            {
                image.Bits[image.Width * yi + xi] = color;
            }
        }

        /// <summary>
        /// Determines if a single point can be drawn on an image (pixel that it would
        /// cover and 4-connected neighbors are equal to the background color).
        /// No subpixel sampling is used.
        /// </summary>
        /// <param name="x">x-coordinate of the point to draw</param>
        /// <param name="y">y-coordinate of the point to draw</param>
        /// <param name="backgColor">Expected background color</param>
        /// <param name="image">Image to draw the point onto</param>
        /// <returns>True if the point can be drawn without rendering over something else or being outside of the image</returns>
        public static bool IsPtFree(double x, double y, int backgColor, DirectBitmap image)
        {
            bool isFree = false;
            int xi = (int)(x + 0.5);
            int yi = (int)(y + 0.5);
            if (xi >= 1 && xi < image.Width - 1 && yi >= 1 && yi < image.Height - 1)
            {
                if (backgColor == image.Bits[image.Width * (yi - 1) + xi] &&
                    backgColor == image.Bits[image.Width * (yi + 1) + xi] &&
                    backgColor == image.Bits[image.Width * yi + xi] &&
                    backgColor == image.Bits[image.Width * yi + xi + 1] &&
                    backgColor == image.Bits[image.Width * yi + xi - 1])
                    isFree = true;
            }
            return isFree;
        }

        /// <summary>
        /// Draws a line on an image. No subpixel sampling is used.
        /// </summary>
        /// <param name="x1">x-coordinate of the first endpoint of the line</param>
        /// <param name="y1">y-coordinate of the first endpoint of the line</param>
        /// <param name="x2">x-coordinate of the second endpoint of the line</param>
        /// <param name="y2">y-coordinate of the second endpoint of the line</param>
        /// <param name="color">Color to give the line</param>
        /// <param name="image">Image to draw the line onto</param>
        public static void DrawLine(double x1, double y1, double x2, double y2, int color, DirectBitmap image)
        {
            int x1i = (int)(x1 + 0.5);
            int y1i = (int)(y1 + 0.5);
            int x2i = (int)(x2 + 0.5);
            int y2i = (int)(y2 + 0.5);
            if (Math.Abs(x2 - x1) > Math.Abs(y2 - y1))
            {
                int xStart = Math.Max(0, Math.Min(x1i, x2i));
                int xEnd = Math.Min(image.Width - 1, Math.Max(x1i, x2i));
                for (int xi = xStart; xi <= xEnd; ++xi)
                    DrawPt(xi, y1 + (xi - x1) * (y2 - y1) / (x2 - x1), color, image);
            }
            else
            {
                int yStart = Math.Max(0, Math.Min(y1i, y2i));
                int yEnd = Math.Min(image.Height - 1, Math.Max(y1i, y2i));
                for (int yi = yStart; yi <= yEnd; ++yi)
                    DrawPt(x1 + (yi - y1) * (x2 - x1) / (y2 - y1), yi, color, image);
            }
        }

        /// <summary>
        /// Determines if a line can be drawn on an image (all pixels that it would
        /// cover are equal to the background color). No subpixel sampling is used.
        /// </summary>
        /// <param name="x1">x-coordinate of the first endpoint of the line</param>
        /// <param name="y1">y-coordinate of the first endpoint of the line</param>
        /// <param name="x2">x-coordinate of the second endpoint of the line</param>
        /// <param name="y2">y-coordinate of the second endpoint of the line</param>
        /// <param name="backgColor">Expected background color</param>
        /// <param name="image">Image where the line is intended to be drawn onto</param>
        /// <returns>True if the line can be drawn without rendering over something else or being outside of the image</returns>
        public static bool IsLineFree(double x1, double y1, double x2, double y2, int backgColor, DirectBitmap image)
        {
            int x1i = (int)(x1 + 0.5);
            int y1i = (int)(y1 + 0.5);
            int x2i = (int)(x2 + 0.5);
            int y2i = (int)(y2 + 0.5);
            if (Math.Abs(x2 - x1) > Math.Abs(y2 - y1))
            {
                int xStart = Math.Max(0, Math.Min(x1i, x2i));
                int xEnd = Math.Min(image.Width - 1, Math.Max(x1i, x2i));
                for (int xi = xStart; xi <= xEnd; ++xi)
                    if (!IsPtFree(xi, y1 + (xi - x1) * (y2 - y1) / (x2 - x1), backgColor, image))
                        return false;
            }
            else
            {
                int yStart = Math.Max(0, Math.Min(y1i, y2i));
                int yEnd = Math.Min(image.Height - 1, Math.Max(y1i, y2i));
                for (int yi = yStart; yi <= yEnd; ++yi)
                    if (!IsPtFree(x1 + (yi - y1) * (x2 - x1) / (y2 - y1), yi, backgColor, image))
                        return false;
            }
            return true;
        }

        /// <summary>
        /// Draws a spiral on an image. No subpixel sampling is used.
        /// </summary>
        /// <param name="x">x-coordiante of the center of the spiral</param>
        /// <param name="y">y-coordinate of the center of the spiral</param>
        /// <param name="thetaInit">Inital theta of the spiral</param>
        /// <param name="rInit">Intial radius of the spiral</param>
        /// <param name="thetaPlusSpan">Span of the spiral in theta</param>
        /// <param name="thetaMinusSpan">Span of the spiral backwards in theta</param>
        /// <param name="thetaStep">Step size in theta</param>
        /// <param name="scalePerRev">Amount the spiral will scale for each revolution</param>
        /// <param name="color">Color to give the spiral</param>
        /// <param name="image">Image to draw the spiral onto</param>
        public static void DrawSpiral(double x, double y, double thetaInit, double rInit, double thetaPlusSpan, double thetaMinusSpan, double thetaStep, double scalePerRev, int color, DirectBitmap image)
        {
            double prevX = double.NaN;
            double prevY = double.NaN;
            for (double baseTheta = -thetaMinusSpan; baseTheta < thetaPlusSpan; baseTheta += thetaStep)
            {
                double theta = thetaInit + baseTheta;
                double revs = baseTheta / (2.0 * Math.PI);
                double alpha = Math.Log(scalePerRev);
                double r = rInit * Math.Exp(alpha * revs);
                double curX = x + r * Math.Cos(theta);
                double curY = y + r * Math.Sin(theta);
                if (!double.IsNaN(prevX) && !double.IsNaN(prevY))
                {
                    DrawLine(prevX, prevY, curX, curY, color, image);
                }
                prevX = curX;
                prevY = curY;
            }
        }

        /// <summary>
        /// Draws a thick spiral on an image. No subpixel sampling is used.
        /// </summary>
        /// <param name="x">x-coordiante of the center of the spiral</param>
        /// <param name="y">y-coordinate of the center of the spiral</param>
        /// <param name="thetaInit">Inital theta of the spiral</param>
        /// <param name="rInit">Intial radius of the spiral</param>
        /// <param name="thetaPlusSpan">Span of the spiral in theta</param>
        /// <param name="thetaMinusSpan">Span of the spiral backwards in theta</param>
        /// <param name="thetaStep">Step size in theta</param>
        /// <param name="scalePerRev">Amount the spiral will scale for each revolution</param>
        /// <param name="relThickness">Relative thickness (ratio of the radius at the current point)</param>
        /// <param name="color">Color to give the spiral</param>
        /// <param name="image">Image to draw the spiral onto</param>
        public static void DrawThickSpiral(double x, double y, double thetaInit, double rInit, double thetaPlusSpan, double thetaMinusSpan, double thetaStep, double scalePerRev, double relThickness, int color, DirectBitmap image)
        {
            double prevX = double.NaN;
            double prevY = double.NaN;
            for (double baseTheta = -thetaMinusSpan; baseTheta < thetaPlusSpan; baseTheta += thetaStep)
            {
                double theta = thetaInit + baseTheta;
                double revs = baseTheta / (2.0 * Math.PI);
                double alpha = Math.Log(scalePerRev);
                double r = rInit * Math.Exp(alpha * revs);
                double curX = x + r * Math.Cos(theta);
                double curY = y + r * Math.Sin(theta);
                if (!double.IsNaN(prevX) && !double.IsNaN(prevY))
                {
                    double halfWidth = 0.5 * relThickness * r;
                    double dX = r * (alpha * Math.Cos(theta) - 2.0 * Math.PI * Math.Sin(theta));
                    double dY = r * (alpha * Math.Sin(theta) + 2.0 * Math.PI * Math.Cos(theta));
                    double mag = Math.Sqrt(dX * dX + dY * dY);
                    dX /= mag;
                    dY /= mag;
                    int x1i = (int)(prevX + 0.5);
                    int y1i = (int)(prevY + 0.5);
                    int x2i = (int)(curX + 0.5);
                    int y2i = (int)(curY + 0.5);
                    if (Math.Abs(curX - prevX) > Math.Abs(curY - prevY))
                    {
                        int xStart = Math.Max(0, Math.Min(x1i, x2i));
                        int xEnd = Math.Min(image.Width - 1, Math.Max(x1i, x2i));
                        for (int xi = xStart; xi <= xEnd; ++xi)
                        {
                            DrawLine(xi + halfWidth * dY, prevY + (xi - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, xi - halfWidth * dY, prevY + (xi - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                            if (xi < xEnd)
                            {
                                DrawLine((xi + 1.0 / 3.0) + halfWidth * dY, prevY + ((xi + 1.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, (xi + 1.0 / 3.0) - halfWidth * dY, prevY + ((xi + 1.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                                DrawLine((xi + 2.0 / 3.0) + halfWidth * dY, prevY + ((xi + 2.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, (xi + 2.0 / 3.0) - halfWidth * dY, prevY + ((xi + 2.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                            }
                        }
                    }
                    else
                    {
                        int yStart = Math.Max(0, Math.Min(y1i, y2i));
                        int yEnd = Math.Min(image.Height - 1, Math.Max(y1i, y2i));
                        for (int yi = yStart; yi <= yEnd; ++yi)
                        {
                            DrawLine(prevX + (yi - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, yi - halfWidth * dX, prevX + (yi - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, yi + halfWidth * dX, color, image);
                            if (yi < yEnd)
                            {
                                DrawLine(prevX + ((yi + 1.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, (yi + 1.0 / 3.0) - halfWidth * dX, prevX + ((yi + 1.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, (yi + 1.0 / 3.0) + halfWidth * dX, color, image);
                                DrawLine(prevX + ((yi + 2.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, (yi + 2.0 / 3.0) - halfWidth * dX, prevX + ((yi + 2.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, (yi + 2.0 / 3.0) + halfWidth * dX, color, image);
                            }
                        }
                    }
                }
                prevX = curX;
                prevY = curY;
            }
        }

        /// <summary>
        /// Draws a tapered thick spiral on an image,
        /// where the thickness shrinks to zero at the endpoints.
        /// No subpixel sampling is used.
        /// </summary>
        /// <param name="x">x-coordiante of the center of the spiral</param>
        /// <param name="y">y-coordinate of the center of the spiral</param>
        /// <param name="thetaInit">Inital theta of the spiral</param>
        /// <param name="rInit">Intial radius of the spiral</param>
        /// <param name="thetaPlusSpan">Span of the spiral in theta</param>
        /// <param name="thetaMinusSpan">Span of the spiral backwards in theta</param>
        /// <param name="thetaStep">Step size in theta</param>
        /// <param name="scalePerRev">Amount the spiral will scale for each revolution</param>
        /// <param name="relThickness">Relative thickness (ratio of the radius at the current point)</param>
        /// <param name="decayParameter">Controls how rapidly the thickness decays towards the endpoints</param>
        /// <param name="color">Color to give the spiral</param>
        /// <param name="image">Image to draw the spiral onto</param>
        public static void DrawTaperedSpiral(double x, double y, double thetaInit, double rInit, double thetaPlusSpan, double thetaMinusSpan, double thetaStep, double scalePerRev, double relThickness, double decayParameter, int color, DirectBitmap image)
        {
            double prevX = double.NaN;
            double prevY = double.NaN;
            for (double baseTheta = -thetaMinusSpan; baseTheta < thetaPlusSpan; baseTheta += thetaStep)
            {
                double theta = thetaInit + baseTheta;
                double revs = baseTheta / (2.0 * Math.PI);
                double alpha = Math.Log(scalePerRev);
                double r = rInit * Math.Exp(alpha * revs);
                double curX = x + r * Math.Cos(theta);
                double curY = y + r * Math.Sin(theta);
                if (!double.IsNaN(prevX) && !double.IsNaN(prevY))
                {
                    double halfWidth = 0.5 * relThickness * r * Math.Exp(-Math.Pow(Math.Tan(-0.5 * Math.PI + Math.PI * (baseTheta + thetaMinusSpan) / (thetaPlusSpan + thetaMinusSpan)) * decayParameter, 2));
                    double dX = r * (alpha * Math.Cos(theta) - 2.0 * Math.PI * Math.Sin(theta));
                    double dY = r * (alpha * Math.Sin(theta) + 2.0 * Math.PI * Math.Cos(theta));
                    double mag = Math.Sqrt(dX * dX + dY * dY);
                    dX /= mag;
                    dY /= mag;
                    int x1i = (int)(prevX + 0.5);
                    int y1i = (int)(prevY + 0.5);
                    int x2i = (int)(curX + 0.5);
                    int y2i = (int)(curY + 0.5);
                    if (Math.Abs(curX - prevX) > Math.Abs(curY - prevY))
                    {
                        int xStart = Math.Max(0, Math.Min(x1i, x2i));
                        int xEnd = Math.Min(image.Width - 1, Math.Max(x1i, x2i));
                        for (int xi = xStart; xi <= xEnd; ++xi)
                        {
                            DrawLine(xi + halfWidth * dY, prevY + (xi - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, xi - halfWidth * dY, prevY + (xi - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                            if (xi < xEnd)
                            {
                                DrawLine((xi + 1.0 / 3.0) + halfWidth * dY, prevY + ((xi + 1.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, (xi + 1.0 / 3.0) - halfWidth * dY, prevY + ((xi + 1.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                                DrawLine((xi + 2.0 / 3.0) + halfWidth * dY, prevY + ((xi + 2.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) - halfWidth * dX, (xi + 2.0 / 3.0) - halfWidth * dY, prevY + ((xi + 2.0 / 3.0) - prevX) * (curY - prevY) / (curX - prevX) + halfWidth * dX, color, image);
                            }
                        }
                    }
                    else
                    {
                        int yStart = Math.Max(0, Math.Min(y1i, y2i));
                        int yEnd = Math.Min(image.Height - 1, Math.Max(y1i, y2i));
                        for (int yi = yStart; yi <= yEnd; ++yi)
                        {
                            DrawLine(prevX + (yi - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, yi - halfWidth * dX, prevX + (yi - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, yi + halfWidth * dX, color, image);
                            if (yi < yEnd)
                            {
                                DrawLine(prevX + ((yi + 1.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, (yi + 1.0 / 3.0) - halfWidth * dX, prevX + ((yi + 1.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, (yi + 1.0 / 3.0) + halfWidth * dX, color, image);
                                DrawLine(prevX + ((yi + 2.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) + halfWidth * dY, (yi + 2.0 / 3.0) - halfWidth * dX, prevX + ((yi + 2.0 / 3.0) - prevY) * (curX - prevX) / (curY - prevY) - halfWidth * dY, (yi + 2.0 / 3.0) + halfWidth * dX, color, image);
                            }
                        }
                    }
                }
                prevX = curX;
                prevY = curY;
            }
        }

        /// <summary>
        /// Draws a spiral on an image. No subpixel sampling is used. Extents of the spiral
        /// are determined by continuing the spiral until it hits the boundary of the image
        /// or a nonzero pixel.
        /// </summary>
        /// <param name="x">x-coordiante of the center of the spiral</param>
        /// <param name="y">y-coordinate of the center of the spiral</param>
        /// <param name="thetaInit">Inital theta of the spiral</param>
        /// <param name="rInit">Intial radius of the spiral</param>
        /// <param name="thetaStep">Step size in theta</param>
        /// <param name="scalePerRev">Amount the spiral will scale for each revolution</param>
        /// <param name="color">Color to give the spiral</param>
        /// <param name="image">Image to draw the spiral onto</param>
        /// <param name="thetaPlusSpan">Computed thetaPlusSpan for the spiral being rendered</param>
        /// <param name="thetaMinusSpan">Computed thetaMinusSpan for the spiral being rendered</param>
        public static void DrawSpiral(double x, double y, double thetaInit, double rInit, double thetaStep, double scalePerRev, int color, DirectBitmap image, out double thetaPlusSpan, out double thetaMinusSpan)
        {
            double prevX = double.NaN;
            double prevY = double.NaN;
            // Determine thetaPlusSpan
            thetaPlusSpan = 0;
            bool keepLooping = true;
            double baseTheta = 0;
            while (keepLooping)
            {
                double theta = thetaInit + baseTheta;
                double revs = baseTheta / (2.0 * Math.PI);
                double alpha = Math.Log(scalePerRev);
                double r = rInit * Math.Exp(alpha * revs);
                double curX = x + r * Math.Cos(theta);
                double curY = y + r * Math.Sin(theta);
                if (!double.IsNaN(prevX) && !double.IsNaN(prevY))
                {
                    if (!IsLineFree(prevX, prevY, curX, curY, unchecked((int)0xff000000), image))
                    {
                        keepLooping = false;
                        thetaPlusSpan = baseTheta;
                    }
                }
                prevX = curX;
                prevY = curY;
                baseTheta += thetaStep;
            }
            prevX = double.NaN;
            prevY = double.NaN;
            // Determine thetaMinusSpan
            thetaMinusSpan = 0;
            keepLooping = true;
            baseTheta = 0;
            while (keepLooping)
            {
                double theta = thetaInit + baseTheta;
                double revs = baseTheta / (2.0 * Math.PI);
                double alpha = Math.Log(scalePerRev);
                double r = rInit * Math.Exp(alpha * revs);
                double curX = x + r * Math.Cos(theta);
                double curY = y + r * Math.Sin(theta);
                if (!double.IsNaN(prevX) && !double.IsNaN(prevY))
                {
                    if (!IsLineFree(prevX, prevY, curX, curY, unchecked((int)0xff000000), image))
                    {
                        keepLooping = false;
                        thetaMinusSpan = -baseTheta;
                    }
                }
                prevX = curX;
                prevY = curY;
                baseTheta -= thetaStep;
            }

            if (thetaMinusSpan > 0 || thetaPlusSpan > 0)
                DrawSpiral(x, y, thetaInit, rInit, thetaPlusSpan, thetaMinusSpan, thetaStep, scalePerRev, color, image);
        }

        /// <summary>
        /// Draws a closed contour onto an image
        /// </summary>
        /// <param name="xCoords">Array of x-coordinates of the contour</param>
        /// <param name="yCoords">Array of y-coordinates of the contour</param>
        /// <param name="color">Color of the contour</param>
        /// <param name="image">Image to draw onto</param>
        public static void DrawClosedContour(double[] xCoords, double[] yCoords, int color, DirectBitmap image)
        {
            int N = xCoords.Length;
            if (yCoords.Length != N)
                throw new Exception("DrawClosedContour: Lengths of xCoords and yCoords do not match!");

            for (int i = 0; i < N; ++i)
            {
                int ii = (i + 1) % N;
                double x1 = xCoords[i];
                double y1 = yCoords[i];
                double x2 = xCoords[ii];
                double y2 = yCoords[ii];
                DrawLine(x1, y1, x2, y2, color, image);
            }
        }

        #endregion Drawing

        #region Filters

        /// <summary>
        /// Inverts a 1-channel image
        /// </summary>
        /// <param name="image">Image to invert</param>
        /// <param name="fixedPt">Fixed point to invert about</param>
        /// <returns>Inverted image</returns>
        public static double[,] Invert(double[,] image, double fixedPt)
        {
            // initialize some constants
            int imW = image.GetLength(1);
            int imH = image.GetLength(0);

            double[,] retVal = new double[imH, imW];

            // Perform the inversion
            for (int y = 0; y < imH; ++y)
            {
                for (int x = 0; x < imW; ++x)
                {
                    retVal[y, x] = fixedPt - image[y, x];
                }
            }

            return retVal;
        }

        /// <summary>
        /// Applies a simple median filter to a 3-channel image of double-precision numbers.
        /// Note that pixels that are too close to the boundary for the full median filter
        /// to operate are left unchanged.
        /// </summary>
        /// <param name="image">Input image</param>
        /// <param name="w">Width of the median filter (must be odd)</param>
        /// <param name="h">Height of the median filter (must be odd)</param>
        /// <returns>Median filtered image</returns>
        public static double[,,] MedianFilter(double[,,] image, int w, int h)
        {
            if (w % 2 != 1 || h % 2 != 1)
                throw new Exception("Width and height of median filter must be odd.");

            // initialize some constants
            int imW = image.GetLength(1);
            int imH = image.GetLength(0);
            int hw = (w - 1) / 2;
            int hh = (h - 1) / 2;

            double[,,] retVal = new double[imH, imW, 3];

            // Perform the filtering
            for (int y = hh; y < imH - hh; ++y)
            {
                for (int x = hw; x < imW - hw; ++x)
                {
                    var channel0Vals = new List<double>();
                    var channel1Vals = new List<double>();
                    var channel2Vals = new List<double>();
                    for (int yofs = -hh; yofs <= hh; ++yofs)
                    {
                        for (int xofs = -hw; xofs <= hw; ++xofs)
                        {
                            channel0Vals.Add(image[y + yofs, x + xofs, 0]);
                            channel1Vals.Add(image[y + yofs, x + xofs, 1]);
                            channel2Vals.Add(image[y + yofs, x + xofs, 2]);
                        }
                    }

                    // Compute medians and store into pixel (x, y)
                    channel0Vals.Sort();
                    channel1Vals.Sort();
                    channel2Vals.Sort();
                    retVal[y, x, 0] = channel0Vals[(channel0Vals.Count - 1) / 2];
                    retVal[y, x, 1] = channel1Vals[(channel1Vals.Count - 1) / 2];
                    retVal[y, x, 2] = channel2Vals[(channel2Vals.Count - 1) / 2];
                }
            }

            return retVal;
        }

        /// <summary>
        /// Applies a simple 3x3 median filter to a 1-channel image of double-precision numbers,
        /// using a fast algorithm (see http://ndevilla.free.fr/median/median.pdf)
        /// Note that pixels that are too close to the boundary for the full median filter
        /// to operate are left unchanged.
        /// </summary>
        /// <param name="image">Input image</param>
        /// <returns>Median filtered image</returns>
        public static double[,] Fast3x3MedianFilter(double[,] image)
        {
            // initialize some constants
            int imW = image.GetLength(1);
            int imH = image.GetLength(0);

            double[,] retVal = new double[imH, imW];

            // Perform the filtering
            int ymax = imH - 1;
            int xmax = imW - 1;
            double v0, v1, v2, v3, v4, v5, v6, v7, v8, tmp;
            for (int y = 1; y < ymax; ++y)
            {
                for (int x = 1; x < xmax; ++x)
                {
                    v0 = image[y - 1, x - 1];
                    v1 = image[y - 1, x];
                    v2 = image[y - 1, x + 1];
                    v3 = image[y, x - 1];
                    v4 = image[y, x];
                    v5 = image[y, x + 1];
                    v6 = image[y + 1, x - 1];
                    v7 = image[y + 1, x];
                    v8 = image[y + 1, x + 1];
                    if (v1 > v2) { tmp = v1; v1 = v2; v2 = tmp; }
                    if (v4 > v5) { tmp = v4; v4 = v5; v5 = tmp; }
                    if (v7 > v8) { tmp = v7; v7 = v8; v8 = tmp; }
                    if (v0 > v1) { tmp = v0; v0 = v1; v1 = tmp; }
                    if (v3 > v4) { tmp = v3; v3 = v4; v4 = tmp; }
                    if (v6 > v7) { tmp = v6; v6 = v7; v7 = tmp; }
                    if (v1 > v2) { tmp = v1; v1 = v2; v2 = tmp; }
                    if (v4 > v5) { tmp = v4; v4 = v5; v5 = tmp; }
                    if (v7 > v8) { tmp = v7; v7 = v8; v8 = tmp; }
                    if (v0 > v3) { tmp = v0; v0 = v3; v3 = tmp; }
                    if (v5 > v8) { tmp = v5; v5 = v8; v8 = tmp; }
                    if (v4 > v7) { tmp = v4; v4 = v7; v7 = tmp; }
                    if (v3 > v6) { tmp = v3; v3 = v6; v6 = tmp; }
                    if (v1 > v4) { tmp = v1; v1 = v4; v4 = tmp; }
                    if (v2 > v5) { tmp = v2; v2 = v5; v5 = tmp; }
                    if (v4 > v7) { tmp = v4; v4 = v7; v7 = tmp; }
                    if (v4 > v2) { tmp = v4; v4 = v2; v2 = tmp; }
                    if (v6 > v4) { tmp = v6; v6 = v4; v4 = tmp; }
                    if (v4 > v2) { tmp = v4; v4 = v2; v2 = tmp; }
                    retVal[y, x] = v4;
                }
            }

            return retVal;
        }

        /// <summary>
        /// Performs a 1-dimensional convolution on 1-d input data, returning the result in the output.
        /// </summary>
        /// <param name="input">Input array</param>
        /// <param name="kernel">Convolution kernel</param>
        /// <param name="shape">Controls the size of the output returned. Options are:<br/>
        /// 0: 'valid' - Return only parts of the convolution that are computed without requiring padding at the boundary.<br/>
        /// 1: 'same' - Return the central part of the convolution, which is the same size as the input.<br/>
        /// 2: 'full' - Return the full 1-d convolution
        /// </param>
        /// <param name="origin">
        /// Element of the kernel to be considered the "origin" (applies when shape == 1).
        /// If -1 is passed in (default value) this number is automatically calculated as (kernel.Length-1)/2
        /// </param>
        /// <param name="boundaryCondition">Specifies how the boundary pixels will be padded (applies when shape > 0). Options are:<br/>
        /// 0: Default value. Pad with the constant value specified in boundaryValue (default = 0.0)
        /// 1: Periodic boundary conditions
        /// 2: Reflect the data at the boundary
        /// </param>
        /// <param name="boundaryValue">Value to use for padding the boundary (if boundaryCondition == 0)</param>
        /// <returns>Result of the convolution</returns>
        public static double[] Convolve1D(double[] input, double[] kernel, int shape, int origin = -1, int boundaryCondition = 0, double boundaryValue = 0.0)
        {
            // initialize some constants
            int imH = input.Length;

            double[] origInput = null;
            int offset = 0; // offset to origInput within new input
            switch (shape)
            {
                case 1:
                    // Pad input sufficiently for same-sized output
                    origInput = input;
                    input = new double[imH + kernel.Length - 1];
                    if (-1 == origin)
                        origin = (kernel.Length - 1) / 2;
                    offset = kernel.Length - 1 - origin;
                    break;
                case 2:
                    // Pad input sufficiently for full-sized output
                    origInput = input;
                    input = new double[imH + 2 * (kernel.Length - 1)];
                    offset = kernel.Length - 1;
                    break;
                default: // case 0, nothing to do
                    break;
            }
            int offsetY = 0;
            offsetY = offset;

            if (null != origInput)
            {
                int hoi = origInput.Length;
                int hi = input.Length;
                int nk = kernel.Length;
                // Copy original input into appropriate section of new input
                for (int y = 0; y < hoi; ++y)
                    input[offsetY + y] = origInput[y];
                // Fill in the padded elements
                switch (boundaryCondition)
                {
                    case 0: // Pad with the constant value specified in boundaryValue (default = 0.0)
                        for (int y = 0; y < offset; ++y)
                            input[y] = boundaryValue;
                        for (int y = offset + hoi; y < hi; ++y)
                            input[y] = boundaryValue;
                        break;
                    case 1: // Periodic boundary conditions
                        for (int y = 0; y < offset; ++y)
                            input[y] = origInput[(y - offset + hoi) % hoi];
                        for (int y = offset + hoi; y < hi; ++y)
                            input[y] = origInput[(y - offset) % hoi];
                        break;
                    case 2: // Reflect the data at the boundary
                        for (int y = 0; y < offset; ++y)
                            input[y] = origInput[offset - 1 - y];
                        for (int y = offset + hoi; y < hi; ++y)
                            input[y] = origInput[2 * (offset + hoi) - 1 - y];
                        break;
                    default:
                        break;
                }
            }

            return Convolve1D(input, kernel);
        }

        /// <summary>
        /// Performs a 1-dimensional convolution on 2-d input data, returning the result in the output.
        /// </summary>
        /// <param name="input">Input array</param>
        /// <param name="kernel">Convolution kernel</param>
        /// <param name="dim">Dimension to convolve along (0 or 1)</param>
        /// <param name="shape">Controls the size of the output returned. Options are:<br/>
        /// 0: 'valid' - Return only parts of the convolution that are computed without requiring padding at the boundary.<br/>
        /// 1: 'same' - Return the central part of the convolution, which is the same size as the input.<br/>
        /// 2: 'full' - Return the full 1-d convolution
        /// </param>
        /// <param name="origin">
        /// Element of the kernel to be considered the "origin" (applies when shape == 1).
        /// If -1 is passed in (default value) this number is automatically calculated as (kernel.Length-1)/2
        /// </param>
        /// <param name="boundaryCondition">Specifies how the boundary pixels will be padded (applies when shape > 0). Options are:<br/>
        /// 0: Default value. Pad with the constant value specified in boundaryValue (default = 0.0)
        /// 1: Periodic boundary conditions
        /// 2: Reflect the data at the boundary
        /// </param>
        /// <param name="boundaryValue">Value to use for padding the boundary (if boundaryCondition == 0)</param>
        /// <returns>Result of the convolution</returns>
        public static double[,] Convolve1D(double[,] input, int dim, double[] kernel, int shape, int origin = -1, int boundaryCondition = 0, double boundaryValue = 0.0)
        {
            // initialize some constants
            int imW = input.GetLength(1);
            int imH = input.GetLength(0);

            double[,] origInput = null;
            int offset = 0; // offset to origInput within new input
            switch (shape)
            {
                case 1:
                    // Pad input sufficiently for same-sized output
                    origInput = input;
                    if (dim == 0)
                        input = new double[imH + kernel.Length - 1, imW];
                    else
                        input = new double[imH, imW + kernel.Length - 1];
                    if (-1 == origin)
                        origin = (kernel.Length - 1) / 2;
                    offset = kernel.Length - 1 - origin;
                    break;
                case 2:
                    // Pad input sufficiently for full-sized output
                    origInput = input;
                    if (dim == 0)
                        input = new double[imH + 2 * (kernel.Length - 1), imW];
                    else
                        input = new double[imH, imW + 2 * (kernel.Length - 1)];
                    offset = kernel.Length - 1;
                    break;
                default: // case 0, nothing to do
                    break;
            }
            int offsetX = 0;
            int offsetY = 0;
            if (dim == 0)
                offsetY = offset;
            else
                offsetX = offset;

            if (null != origInput)
            {
                int woi = origInput.GetLength(1);
                int hoi = origInput.GetLength(0);
                int wi = input.GetLength(1);
                int hi = input.GetLength(0);
                int nk = kernel.Length;
                // Copy original input into appropriate section of new input
                for (int y = 0; y < hoi; ++y)
                    for (int x = 0; x < woi; ++x)
                        input[offsetY + y, offsetX + x] = origInput[y, x];
                // Fill in the padded elements
                if (dim == 0)
                {
                    switch (boundaryCondition)
                    {
                        case 0: // Pad with the constant value specified in boundaryValue (default = 0.0)
                            for (int y = 0; y < offset; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = boundaryValue;
                            for (int y = offset + hoi; y < hi; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = boundaryValue;
                            break;
                        case 1: // Periodic boundary conditions
                            for (int y = 0; y < offset; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = origInput[(y - offset + hoi) % hoi, x];
                            for (int y = offset + hoi; y < hi; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = origInput[(y - offset) % hoi, x];
                            break;
                        case 2: // Reflect the data at the boundary
                            for (int y = 0; y < offset; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = origInput[offset - 1 - y, x];
                            for (int y = offset + hoi; y < hi; ++y)
                                for (int x = 0; x < wi; ++x)
                                    input[y, x] = origInput[2 * (offset + hoi) - 1 - y, x];
                            break;
                        default:
                            break;
                    }
                }
                else
                {
                    switch (boundaryCondition)
                    {
                        case 0: // Pad with the constant value specified in boundaryValue (default = 0.0)
                            for (int y = 0; y < hi; ++y)
                                for (int x = 0; x < offset; ++x)
                                    input[y, x] = boundaryValue;
                            for (int y = 0; y < hi; ++y)
                                for (int x = offset + woi; x < wi; ++x)
                                    input[y, x] = boundaryValue;
                            break;
                        case 1: // Periodic boundary conditions
                            for (int y = 0; y < hi; ++y)
                                for (int x = 0; x < offset; ++x)
                                    input[y, x] = origInput[y, (x - offset + woi) % woi];
                            for (int y = 0; y < hi; ++y)
                                for (int x = offset + woi; x < wi; ++x)
                                    input[y, x] = origInput[y, (x - offset) % woi];
                            break;
                        case 2: // Reflect the data at the boundary
                            for (int y = 0; y < hi; ++y)
                                for (int x = 0; x < offset; ++x)
                                    input[y, x] = origInput[y, offset - 1 - x];
                            for (int y = 0; y < hi; ++y)
                                for (int x = offset + woi; x < wi; ++x)
                                    input[y, x] = origInput[y, 2 * (offset + woi) - 1 - x];
                            break;
                        default:
                            break;
                    }
                }
            }

            return Convolve1D(input, dim, kernel);
        }

        /// <summary>
        /// Construct a gaussian kernel normalized to weight 1.
        /// </summary>
        /// <param name="sigma">Standard deviation of the gaussian kernel, in units of the index into the array.</param>
        /// <param name="halfWidth">Desired window half-width. Number of elements in the returned array will equal 2 * halfWidth + 1</param>
        /// <param name="smoothOffset">
        /// If true, subtract a constant from the gaussian function such that the next sample on either side of the window will equal zero.
        /// Since we cannot represent a true gaussian function (which has infinite support), this allows us to generate a finite kernel
        /// that has a smoother transition to zero at the boundary.
        /// </param>
        /// <returns>Gaussian kernel with the desired standard deviation and window half-width. Number of elements will equal 2 * halfWidth + 1</returns>
        public static double[] ConstructGaussianKernel(double sigma, int halfWidth, bool smoothOffset = false)
        {
            double[] output = new double[2 * halfWidth + 1];
            for (int i = 0; i <= halfWidth; ++i)
            {
                output[halfWidth + i] = Math.Exp(-i * i / (2.0 * sigma * sigma));
                output[halfWidth - i] = output[halfWidth + i];
            }
            if (smoothOffset)
            {
                double offset = Math.Exp(-(halfWidth + 1) * (halfWidth + 1) / (2.0 * sigma * sigma));
                for (int i = 0; i < 2 * halfWidth + 1; ++i)
                    output[i] -= offset;
            }
            // normalize
            double sum = 0.0;
            for (int i = 0; i < 2 * halfWidth + 1; ++i)
                sum += output[i];
            for (int i = 0; i < 2 * halfWidth + 1; ++i)
                output[i] /= sum;
            return output;
        }

        #region Private Helper Functions

        /// <summary>
        /// Helper function to perform a 1-dimensional convolution on 1-d input data, returning the result in the output.
        /// This function is equivalent to calling Convolve1D(input, kernel, 0), computing only the 'valid' portion of the
        /// convolution (i.e. those values that are computed without requiring padding at the boundary).<br/>
        /// </summary>
        /// <param name="input">Input array</param>
        /// <param name="kernel">Convolution kernel</param>
        /// <returns>Result of the convolution</returns>
        private static double[] Convolve1D(double[] input, double[] kernel)
        {
            double[] output = new double[0];
            int ni = input.Length;
            int nk = kernel.Length;
            int no = ni - nk + 1;
            if (no <= 0)
                return new double[0];
            output = new double[no];

            for (int i = 0; i < no; ++i)
            {
                for (int j = 0; j < nk; ++j)
                    output[i] += input[i + nk - 1 - j] * kernel[j];
            }
            return output;
        }

        /// <summary>
        /// Helper function to perform a 1-dimensional convolution on 2-d input data, returning the result in the output.
        /// This function is equivalent to calling Convolve1D(input, kernel, 0), computing only the 'valid' portion of the
        /// convolution (i.e. those values that are computed without requiring padding at the boundary).<br/>
        /// </summary>
        /// <param name="input">Input array</param>
        /// <param name="dim">Dimension to convolve along (0 or 1)</param>
        /// <param name="kernel">Convolution kernel</param>
        /// <returns>Result of the convolution</returns>
        private static double[,] Convolve1D(double[,] input, int dim, double[] kernel)
        {
            double[,] output = new double[0, 0];
            if (dim == 0)
            {
                int ni = input.GetLength(0);
                int nk = kernel.Length;
                int no = ni - nk + 1;
                if (no <= 0)
                    return new double[0, 0];
                int w = input.GetLength(1);
                output = new double[no, w];
                for (int x = 0; x < w; ++x)
                {
                    for (int i = 0; i < no; ++i)
                    {
                        for (int j = 0; j < nk; ++j)
                            output[i, x] += input[i + nk - 1 - j, x] * kernel[j];
                    }
                }
            }
            else
            {
                int ni = input.GetLength(1);
                int nk = kernel.Length;
                int no = ni - nk + 1;
                if (no <= 0)
                    return new double[0, 0];
                int h = input.GetLength(0);
                output = new double[h, no];
                for (int y = 0; y < h; ++y)
                {
                    for (int i = 0; i < no; ++i)
                    {
                        for (int j = 0; j < nk; ++j)
                            output[y, i] += input[y, i + nk - 1 - j] * kernel[j];
                    }
                }
            }
            return output;
        }

        #endregion Private Helper Functions

        #endregion Filters

        #region Segmentation

        /// <summary>
        /// Simple thresholding operation
        /// </summary>
        /// <param name="inputImage">Input image</param>
        /// <param name="threshold">Threshold to apply</param>
        /// <param name="fillValue">Fill value to apply for points that pass the threshold (default = 255)</param>
        /// <param name="type">Type of threshold (default = 0). Four types are supported:<br/>
        /// 0: pixel &gt; threshold<br/>
        /// 1: pixel &lt; threshold<br/>
        /// 2: pixel &gt;= threshold<br/>
        /// 3: pixel &lt;= threshold<br/>
        /// </param>
        /// <returns>Image where fillValue denotes values that pass the threshold and 0 denotes values that fail</returns>
        public static byte[,] Threshold(double[,] inputImage, double threshold, byte fillValue = 255, int type = 0)
        {
            // initialization
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            byte[,] outputImage = new byte[h, w];
            switch (type)
            {
                case 0:
                    for (int y = 0; y < h; ++y)
                        for (int x = 0; x < w; ++x)
                            if (inputImage[y, x] > threshold)
                                outputImage[y, x] = fillValue;
                    break;
                case 1:
                    for (int y = 0; y < h; ++y)
                        for (int x = 0; x < w; ++x)
                            if (inputImage[y, x] < threshold)
                                outputImage[y, x] = fillValue;
                    break;
                case 2:
                    for (int y = 0; y < h; ++y)
                        for (int x = 0; x < w; ++x)
                            if (inputImage[y, x] >= threshold)
                                outputImage[y, x] = fillValue;
                    break;
                case 3:
                    for (int y = 0; y < h; ++y)
                        for (int x = 0; x < w; ++x)
                            if (inputImage[y, x] <= threshold)
                                outputImage[y, x] = fillValue;
                    break;
                default:
                    break;
            }
            return outputImage;
        }

        /// <summary>
        /// Masks an image with a binary mask (nonzero pixels are interpreted as "true")
        /// </summary>
        /// <param name="image">Original image</param>
        /// <param name="mask">Mask to apply</param>
        /// <returns>Masked image</returns>
        public static DirectBitmap Mask(DirectBitmap image, byte[,] mask)
        {
            int w = image.Width;
            int h = image.Height;
            DirectBitmap retVal = new DirectBitmap(w, h);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    retVal.Bits[w * y + x] = mask[y, x] != 0 ? image.Bits[w * y + x] : unchecked((int)0xff000000);
            return retVal;
        }

        /// <summary>
        /// Masks an image with a binary mask (nonzero pixels are interpreted as "true" if invert is not enabled,
        /// and "false" if invert is enabled).
        /// </summary>
        /// <param name="image">Original image</param>
        /// <param name="mask">Mask to apply</param>
        /// <param name="invert">If true, invert the mask. If invert is true, pixels equal to 0 are considered part of the mask and nonzero pixels are not.</param>
        /// <param name="background">Background color to use for the pixels that are not part of the mask</param>
        /// <returns>Masked image</returns>
        public static DirectBitmap Mask(DirectBitmap image, byte[,] mask, bool invert, int background)
        {
            int w = image.Width;
            int h = image.Height;
            DirectBitmap retVal = new DirectBitmap(w, h);
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    retVal.Bits[w * y + x] = (mask[y, x] != 0 ^ invert) ? image.Bits[w * y + x] : background;
            return retVal;
        }

        /// <summary>
        /// Masks an image with a binary mask (nonzero pixels are interpreted as "true")
        /// </summary>
        /// <param name="image">Original image</param>
        /// <param name="mask">Mask to apply</param>
        /// <returns>Masked image</returns>
        public static byte[,] Mask(byte[,] image, byte[,] mask)
        {
            int w = image.GetLength(1);
            int h = image.GetLength(0);
            byte[,] retVal = new byte[h, w];
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    retVal[y, x] = mask[y, x] != 0 ? image[y, x] : (byte)0;
            return retVal;
        }

        #endregion Segmentation

        #region Binary Morph Ops

        /// <summary>
        /// Dilates a binary image by a given kernel (pixels are assumed to be either 0 or 255).
        /// Pixels outside of the boundary are assumed to equal zero.
        /// </summary>
        /// <param name="inputImage">Input image</param>
        /// <param name="kernel">Dilation kernel</param>
        /// <param name="rowOrigin">Row that should be considered the origin in the kernel</param>
        /// <param name="colOrigin">Column that should be considered the origin in the kernel</param>
        /// <returns>Dilated image</returns>
        public static byte[,] BWDilate(byte[,] inputImage, bool[,] kernel, int rowOrigin, int colOrigin)
        {
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            byte[,] result = new byte[h, w];
            int wk = kernel.GetLength(1);
            int hk = kernel.GetLength(0);

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    bool found = false;
                    for (int yofs = Math.Max(-rowOrigin, -y); yofs < Math.Min(hk - rowOrigin, h - y); ++yofs)
                    {
                        for (int xofs = Math.Max(-colOrigin, -x); xofs < Math.Min(wk - colOrigin, w - x); ++xofs)
                        {
                            if (kernel[rowOrigin + yofs, colOrigin + xofs] && inputImage[y + yofs, x + xofs] != 0)
                            {
                                found = true;
                                break;
                            }
                        }
                        if (found)
                            break;
                    }
                    if (found)
                        result[y, x] = 255;
                }
            }

            return result;
        }

        /// <summary>
        /// Erodes a binary image by a given kernel (pixels are assumed to be either 0 or 255).
        /// Pixels outside of the boundary are assumed to equal 255.
        /// </summary>
        /// <param name="inputImage">Input image</param>
        /// <param name="kernel">Dilation kernel</param>
        /// <param name="rowOrigin">Row that should be considered the origin in the kernel</param>
        /// <param name="colOrigin">Column that should be considered the origin in the kernel</param>
        /// <returns>Eroded image</returns>
        public static byte[,] BWErode(byte[,] inputImage, bool[,] kernel, int rowOrigin, int colOrigin)
        {
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            byte[,] result = new byte[h, w];
            int wk = kernel.GetLength(1);
            int hk = kernel.GetLength(0);

            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    bool found = false;
                    for (int yofs = Math.Max(-rowOrigin, -y); yofs < Math.Min(hk - rowOrigin, h - y); ++yofs)
                    {
                        for (int xofs = Math.Max(-colOrigin, -x); xofs < Math.Min(wk - colOrigin, w - x); ++xofs)
                        {
                            if (kernel[rowOrigin + yofs, colOrigin + xofs] && inputImage[y + yofs, x + xofs] == 0)
                            {
                                found = true;
                                break;
                            }
                        }
                        if (found)
                            break;
                    }
                    if (!found)
                        result[y, x] = 255;
                }
            }

            return result;
        }

        /// <summary>
        /// Subtracts one binary image from another
        /// </summary>
        /// <param name="inputImage">Image to subtract from</param>
        /// <param name="imageToSubtract">Image to subtract</param>
        /// <returns>Result of the subtraction</returns>
        public static byte[,] BWSubtract(byte[,] inputImage, byte[,] imageToSubtract)
        {
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            byte[,] result = new byte[h, w];
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    result[y, x] = (inputImage[y, x] != 0 && imageToSubtract[y, x] == 0) ? inputImage[y, x] : (byte)0;
            return result;
        }

        /// <summary>
        ///     Mimic the behavior of Matlab's bwselect function, except that the input
        ///     image is not binary. Instead, the algorithm simulates floodfilling the region
        ///     equal to the color specified by the user, returning the region that would
        ///     have been filled in the output image (backgColor = not filled, foregColor = filled)
        /// </summary>
        /// <param name="inputImage">Input image</param>
        /// <param name="outputImage">Output image</param>
        /// <param name="seeds">List of seed points (x, y) at which to begin the floodfilling operation</param>
        /// <param name="oldColor">grayscale value of unlabeled blob</param>
        /// <param name="backgColor">grayscale value to use as the background of the output image</param>
        /// <param name="foregColor">grayscale value to use as the foreground of the output image</param>
        /// <param name="contour">Stores the (x, y) positions that were filled</param>
        /// <remarks>
        /// Algorithm uses a queue data structure to avoid having to use recursion.
        /// </remarks>
        public static void BWSelect(byte[,] inputImage, out byte[,] outputImage, List<Tuple<int, int>> seeds, byte oldColor, byte backgColor, byte foregColor, out List<Tuple<int, int>> contour)
        {
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            byte[,] result = new byte[h, w];

            int qStart = 0;
            int x, y;

            // Initialize the output image
            outputImage = new byte[h, w];
            for (y = 0; y < h; ++y)
                for (x = 0; x < w; ++x)
                    outputImage[y, x] = backgColor;

            // Initialize the contour
            contour = new List<Tuple<int, int>>();

            if (backgColor == foregColor)
                return; // nothing to do

            // Initialize the queue, filling the pixels with the 
            // new color immediately so they can never again enter the queue
            for (int idx = 0; idx < seeds.Count; ++idx)
            {
                int x0 = seeds[idx].Item1;
                int y0 = seeds[idx].Item2;
                if (inputImage[y0, x0] == oldColor && outputImage[y0, x0] != foregColor)
                {
                    // push (x,y) onto the queue
                    contour.Add(Tuple.Create(x0, y0));
                    // fill pixel (x0,y0)
                    outputImage[y0, x0] = foregColor;
                }
            }

            while (qStart < contour.Count)
            {
                // pop (x,y) off of the queue
                x = contour[qStart].Item1;
                y = contour[qStart].Item2;
                qStart += 1;

                // push any of the 8-connected neighbors onto the queue, if their color matches oldColor
                // and in the new image they are not already set to foregColor,
                // filling the pixel with the new color immediately so it can never again enter the queue
                if (y > 0)
                {
                    if (x > 0 && inputImage[y - 1, x - 1] == oldColor && outputImage[y - 1, x - 1] != foregColor)
                    {
                        // push (x-1,y-1) onto the queue
                        contour.Add(Tuple.Create(x - 1, y - 1));
                        // fill pixel (x-1,y-1)
                        outputImage[y - 1, x - 1] = foregColor;
                    }
                    if (inputImage[y - 1, x] == oldColor && outputImage[y - 1, x] != foregColor)
                    {
                        // push (x,y-1) onto the queue
                        contour.Add(Tuple.Create(x, y - 1));
                        // fill pixel (x,y-1)
                        outputImage[y - 1, x] = foregColor;
                    }
                    if (x < w - 1 && inputImage[y - 1, x + 1] == oldColor && outputImage[y - 1, x + 1] != foregColor)
                    {
                        // push (x+1,y-1) onto the queue
                        contour.Add(Tuple.Create(x + 1, y - 1));
                        // fill pixel (x+1,y-1)
                        outputImage[y - 1, x + 1] = foregColor;
                    }
                }
                if (x > 0 && inputImage[y, x - 1] == oldColor && outputImage[y, x - 1] != foregColor)
                {
                    // push (x-1,y) onto the queue
                    contour.Add(Tuple.Create(x - 1, y));
                    // fill pixel (x-1,y)
                    outputImage[y, x - 1] = foregColor;
                }
                if (x < w - 1 && inputImage[y, x + 1] == oldColor && outputImage[y, x + 1] != foregColor)
                {
                    // push (x+1,y) onto the queue
                    contour.Add(Tuple.Create(x + 1, y));
                    // fill pixel (x+1,y)
                    outputImage[y, x + 1] = foregColor;
                }
                if (y < h - 1)
                {
                    if (x > 0 && inputImage[y + 1, x - 1] == oldColor && outputImage[y + 1, x - 1] != foregColor)
                    {
                        // push (x-1,y+1) onto the queue
                        contour.Add(Tuple.Create(x - 1, y + 1));
                        // fill pixel (x-1,y+1)
                        outputImage[y + 1, x - 1] = foregColor;
                    }
                    if (inputImage[y + 1, x] == oldColor && outputImage[y + 1, x] != foregColor)
                    {
                        // push (x,y+1) onto the queue
                        contour.Add(Tuple.Create(x, y + 1));
                        // fill pixel (x,y+1)
                        outputImage[y + 1, x] = foregColor;
                    }
                    if (x < w - 1 && inputImage[y + 1, x + 1] == oldColor && outputImage[y + 1, x + 1] != foregColor)
                    {
                        // push (x+1,y+1) onto the queue
                        contour.Add(Tuple.Create(x + 1, y + 1));
                        // fill pixel (x+1,y+1)
                        outputImage[y + 1, x + 1] = foregColor;
                    }
                }
            }
        }

        /// <summary>
        /// Simple boundary detection. Render in white all pixels where any
        /// neighbors differed from the pixel at the center.
        /// </summary>
        /// <param name="inputImage">Original image</param>
        /// <param name="destinationImage">Destination image for the boundary (must be the same size as inputImage)</param>
        public static void BWBoundary(DirectBitmap inputImage, DirectBitmap destinationImage)
        {
            int w = inputImage.Width;
            int h = inputImage.Height;
            if (destinationImage.Width != w || destinationImage.Height != h)
                throw new Exception("BWBoundary: size of destination image does not match size of input image!");

            int white = unchecked((int)0xffffffff); // white
            int black = unchecked((int)0xff000000); // black

            // Set image borders to black
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; x += w - 1)
                    destinationImage.Bits[w * y + x] = black;
            for (int y = 0; y < h; y += h - 1)
                for (int x = 0; x < w; ++x)
                    destinationImage.Bits[w * y + x] = black;

            for (int y = 1; y < h - 1; ++y)
            {
                for (int x = 1; x < w - 1; ++x)
                {
                    int cntr = inputImage.Bits[w * y + x];
                    if (inputImage.Bits[w * (y - 1) + x - 1] == cntr &&
                        inputImage.Bits[w * (y - 1) + x] == cntr &&
                        inputImage.Bits[w * (y - 1) + x + 1] == cntr &&
                        inputImage.Bits[w * y + x - 1] == cntr &&
                        inputImage.Bits[w * y + x + 1] == cntr &&
                        inputImage.Bits[w * (y + 1) + x - 1] == cntr &&
                        inputImage.Bits[w * (y + 1) + x] == cntr &&
                        inputImage.Bits[w * (y + 1) + x + 1] == cntr)
                        destinationImage.Bits[w * y + x] = black;
                    else
                        destinationImage.Bits[w * y + x] = white;
                }
            }
        }

        /// <summary>
        /// Simple boundary detection. Render in white all pixels where any
        /// neighbors differed from the pixel at the center.
        /// </summary>
        /// <param name="inputImage">Original image</param>
        /// <param name="destinationImage">Destination image for the boundary (must be the same size as inputImage)</param>
        /// <param name="restrictColor">Only label as boundary pixels, pixels which are this color</param>
        public static void BWBoundary(DirectBitmap inputImage, DirectBitmap destinationImage, int restrictColor)
        {
            int w = inputImage.Width;
            int h = inputImage.Height;
            if (destinationImage.Width != w || destinationImage.Height != h)
                throw new Exception("BWBoundary: size of destination image does not match size of input image!");

            int white = unchecked((int)0xffffffff); // white
            int black = unchecked((int)0xff000000); // black

            // Set image borders to black
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; x += w - 1)
                    destinationImage.Bits[w * y + x] = black;
            for (int y = 0; y < h; y += h - 1)
                for (int x = 0; x < w; ++x)
                    destinationImage.Bits[w * y + x] = black;

            for (int y = 1; y < h - 1; ++y)
            {
                for (int x = 1; x < w - 1; ++x)
                {
                    int cntr = inputImage.Bits[w * y + x];
                    if ((inputImage.Bits[w * (y - 1) + x - 1] == cntr &&
                        inputImage.Bits[w * (y - 1) + x] == cntr &&
                        inputImage.Bits[w * (y - 1) + x + 1] == cntr &&
                        inputImage.Bits[w * y + x - 1] == cntr &&
                        inputImage.Bits[w * y + x + 1] == cntr &&
                        inputImage.Bits[w * (y + 1) + x - 1] == cntr &&
                        inputImage.Bits[w * (y + 1) + x] == cntr &&
                        inputImage.Bits[w * (y + 1) + x + 1] == cntr) ||
                        inputImage.Bits[w * y + x] != restrictColor)
                        destinationImage.Bits[w * y + x] = black;
                    else
                        destinationImage.Bits[w * y + x] = white;
                }
            }
        }

        /// <summary>
        /// Simple boundary detection (4-connected). Render in white all pixels where any
        /// neighbors differed from the pixel at the center.
        /// </summary>
        /// <param name="inputImage">Original image</param>
        /// <param name="destinationImage">Destination image for the boundary (must be the same size as inputImage)</param>
        /// <param name="restrictColor">Only label as boundary pixels, pixels which are this color</param>
        public static void BWBoundary4(byte[,] inputImage, byte[,] destinationImage, byte restrictColor)
        {
            int w = inputImage.GetLength(1);
            int h = inputImage.GetLength(0);
            if (destinationImage.GetLength(1) != w || destinationImage.GetLength(0) != h)
                throw new Exception("BWBoundary: size of destination image does not match size of input image!");

            // Set image borders to black
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; x += w - 1)
                    destinationImage[y, x] = 0;
            for (int y = 0; y < h; y += h - 1)
                for (int x = 0; x < w; ++x)
                    destinationImage[y, x] = 0;

            for (int y = 1; y < h - 1; ++y)
            {
                for (int x = 1; x < w - 1; ++x)
                {
                    int cntr = inputImage[y, x];
                    if ((inputImage[y - 1, x] == cntr &&
                        inputImage[y, x - 1] == cntr &&
                        inputImage[y, x + 1] == cntr &&
                        inputImage[y + 1, x] == cntr) ||
                        inputImage[y, x] != restrictColor)
                        destinationImage[y, x] = 0;
                    else
                        destinationImage[y, x] = 255;
                }
            }
        }

        /// <summary>
        /// Inverts an image
        /// </summary>
        /// <param name="inputImage">Image to invert</param>
        public static void BWInvert(DirectBitmap inputImage)
        {
            int w = inputImage.Width;
            int h = inputImage.Height;
            for (int y = 0; y < h; ++y)
                for (int x = 0; x < w; ++x)
                    inputImage.Bits[w * y + x] = unchecked((int)0xff000000) | ~inputImage.Bits[w * y + x];
        }

        #endregion Binary Morph Ops

        #region Contour Following

        /// <summary>
        /// Extracts a single closed contour from a binary image. It is assumed that
        /// the contour is only 1 pixel wide, so that the ordering of the pixels
        /// around the contour is unambiguous. Note that the image is modified in this
        /// process (contour is replaced with grey).
        /// </summary>
        /// <param name="bwImage">Input image of edge pixels</param>
        /// <param name="xCoords">X-coordinates of the contour</param>
        /// <param name="yCoords">Y-coordinates of the contour</param>
        public static void FindClosedContour(byte[,] bwImage, out double[] xCoords, out double[] yCoords)
        {
            xCoords = null;
            yCoords = null;

            int w = bwImage.GetLength(1);
            int h = bwImage.GetLength(0);

            // Search the image for a nonzero pixel
            int seedX = -1;
            int seedY = -1;
            for (int y = 0; y < h; ++y)
            {
                for (int x = 0; x < w; ++x)
                {
                    if (bwImage[y, x] != 0)
                    {
                        seedX = x;
                        seedY = y;
                        break;
                    }
                }
                if (seedX >= 0 && seedY >= 0)
                    break;
            }

            if (seedX < 0 || seedY < 0)
                return;

            // Traverse the contour, starting from initPt and terminating if either
            // there is no next point or we come back to initPt.
            var curPt = Tuple.Create(seedX, seedY);
            bwImage[curPt.Item2, curPt.Item1] = 128; // Mark the current point as visited
            var xValues = new List<double>();
            var yValues = new List<double>();
            xValues.Add(curPt.Item1);
            yValues.Add(curPt.Item2); // Add to contour
            bool atStart = true;
            
            while (true) // we explicitly break out of the loop when there are no candidate next pixels
            {
                // Create list of candidate next pixels
                var candidateNextPts = new List<Tuple<int, int>>();
                var nxtPt = Tuple.Create(curPt.Item1 - 1, curPt.Item2 - 1);
                if (nxtPt.Item1 >= 0 && nxtPt.Item2 >= 0 && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1, curPt.Item2 - 1);
                if (nxtPt.Item2 >= 0 && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1 + 1, curPt.Item2 - 1);
                if (nxtPt.Item1 < w && nxtPt.Item2 >= 0 && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1 - 1, curPt.Item2);
                if (nxtPt.Item1 >= 0 && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1 + 1, curPt.Item2);
                if (nxtPt.Item1 < w && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1 - 1, curPt.Item2 + 1);
                if (nxtPt.Item1 >= 0 && nxtPt.Item2 < h && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1, curPt.Item2 + 1);
                if (nxtPt.Item2 < h && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);
                nxtPt = Tuple.Create(curPt.Item1 + 1, curPt.Item2 + 1);
                if (nxtPt.Item1 < w && nxtPt.Item2 < h && bwImage[nxtPt.Item2, nxtPt.Item1] == 255)
                    candidateNextPts.Add(nxtPt);

                if (!atStart && candidateNextPts.Count > 1)
                    throw new Exception("FindClosedContour can only work on contours that are 1-pixel wide (with respect to 8-connected neighbors)");

                if (candidateNextPts.Count == 0)
                    break; // exit loop if there are no next points

                atStart = false; // no longer at the start of the contour
                nxtPt = candidateNextPts[0];

                // Set the current point to the next point and mark as visited, adding it to the list
                curPt = nxtPt;
                bwImage[curPt.Item2, curPt.Item1] = 128; // Mark the current point as visited
                xValues.Add(curPt.Item1);
                yValues.Add(curPt.Item2); // Add to contour
            }

            xCoords = xValues.ToArray();
            yCoords = yValues.ToArray();
        }

        /// <summary>
        ///     Apply a curvature-based smoothing algorithm to the contour.
        /// </summary>
        /// <param name="xCoords">Array of x-coordinates of the contour</param>
        /// <param name="yCoords">Array of y-coordinates of the contour</param>
        /// <param name="smoothingTime">Total smoothing "time".</param>
        /// <param name="maxTStep">Minimum time-step to use</param>
        /// <param name="maxTStep">Maximum time-step to use</param>
        /// <param name="isPeriodic">Whether or not to use periodic boundary conditions</param>
        /// <remarks>
        /// Algorithm moves each point in the direction of the normal at a
        /// velocity equal to the curvature (in radians per unit length). The smoothing "speed" can
        /// be controlled by maxTStep and the total amount of smoothing can be controlled by smoothingTime.
        /// </remarks>
        public static void CurvatureSmoothContour(double[] xCoords, double[] yCoords, double smoothingTime, double minTStep, double maxTStep, bool isPeriodic)
        {
            double epsilon = 1.0e-6f; // positive number very close to zero
            double huge = 1.0e9f; // very large positive number

            int n = xCoords.Length;

            if (n <= 2) // Can't smooth a component consisting of just 2 points
                return;

            // Apply the curvature-driven smoothing to the contour
            var tempX = new double[n];
            var tempY = new double[n];
            var curvatures = new double[n];
            double x0, y0, x1, y1, x2, y2, xt, yt, magT, magN;
            double t = 0.0f, tStep;
            bool running = true;
            while (running)
            {
                #region - Find tStep, and compute curvatures and normals -

                tStep = maxTStep;
                for (int p = 0; p < n; p++)
                {
                    x0 = xCoords[(p + n - 1) % n];
                    y0 = yCoords[(p + n - 1) % n];
                    x1 = xCoords[p];
                    y1 = yCoords[p];
                    x2 = xCoords[(p + 1) % n];
                    y2 = yCoords[(p + 1) % n];

                    // Use zero-curvature boundary conditions if periodic is not specified
                    if (!isPeriodic && (p == 0 || p == n - 1))
                    {
                        curvatures[p] = 0.0;
                        tempX[p] = 0.0;
                        tempY[p] = 0.0;
                        continue;
                    }
                    xt = 0.5 * (x2 - x0);
                    yt = 0.5 * (y2 - y0); // (xt, yt) is the approximate tangent vector
                    magT = Math.Sqrt(xt * xt + yt * yt);
                    if (magT > epsilon)
                        curvatures[p] = Math.Abs(ThetaChange(x1 - x0, y1 - y0, x2 - x1, y2 - y1)) / magT;
                    else
                        curvatures[p] = huge;
                    // If the curvature is too small, don't try to move the point
                    if (curvatures[p] <= epsilon)
                    {
                        tempX[p] = 0.0;
                        tempY[p] = 0.0;
                        continue;
                    }
                    tempX[p] = 0.5 * (x2 + x0) - x1;
                    tempY[p] = 0.5 * (y2 + y0) - y1;
                    // Use (tempX[p], tempY[p]) to store the approximate normal vector
                    magN = Math.Sqrt(tempX[p] * tempX[p] + tempY[p] * tempY[p]);
                    if (magN <= epsilon)
                    {
                        // If the normal is unstable, don't try to move the point
                        tempX[p] = 0.0;
                        tempY[p] = 0.0;
                        continue;
                    }

                    // Limit the motion to the size of magN, so the motion of the curve doesn't over-shoot
                    // the average of the neighbors (unless curvture == huge, in which case tStep will be set too small unless we stop it).
                    if (magN < curvatures[p] * tStep && curvatures[p] != huge)
                        tStep = magN / curvatures[p];
                } // end for (int p = 0; p < n; p++)
                if (tStep > maxTStep)
                    tStep = maxTStep;
                if (tStep < minTStep)
                    tStep = minTStep;
                if (smoothingTime < t + tStep)
                {
                    tStep = smoothingTime - t;
                    if (tStep < 0.0f) // In case of floating-point rounding error
                        tStep = 0.0f;
                    running = false;
                }

                #endregion - Find tStep, and compute curvatures and normals -

                #region - Apply smoothing -

                for (int p = 0; p < n; p++)
                {
                    magN = Math.Sqrt(tempX[p] * tempX[p] + tempY[p] * tempY[p]);
                    if (magN > epsilon)
                    {
                        if (tStep * curvatures[p] < magN)
                        // ensure the point doesn't overshoot the midpoint of the neighbors
                        {
                            // Add tStep * curvatures[p] * unit-normal to the position of the point
                            tempX[p] = xCoords[p] + tStep * curvatures[p] * tempX[p] / magN;
                            tempY[p] = yCoords[p] + tStep * curvatures[p] * tempY[p] / magN;
                        }
                        else
                        {
                            // Move the point all the way to the midpoint of the neighbors
                            tempX[p] = xCoords[p] + tempX[p];
                            tempY[p] = yCoords[p] + tempY[p];
                        }
                    }
                    else
                    {
                        // If the normal is unstable, don't try to move the point
                        tempX[p] = xCoords[p];
                        tempY[p] = yCoords[p];
                    }
                }
                for (int p = 0; p < n; p++)
                {
                    xCoords[p] = tempX[p];
                    yCoords[p] = tempY[p];
                }

                #endregion - Apply smoothing -

                t += tStep;
            } // end while (running)
        } // end method CurvatureSmoothContour

        /// <summary>
        ///     Compute the change in theta from one vector to the next
        /// </summary>
        /// <param name="x0">x-component of the first vector</param>
        /// <param name="y0">y-component of the first vector</param>
        /// <param name="x1">x-component of the second vector</param>
        /// <param name="y1">y-component of the second vector</param>
        /// <returns>Change in theta</returns>
        private static double ThetaChange(double x0, double y0, double x1, double y1)
        {
            double epsilon = 1.0e-6;
            double mag0 = Math.Sqrt(x0 * x0 + y0 * y0);
            double mag1 = Math.Sqrt(x1 * x1 + y1 * y1);
            if (mag0 * mag1 <= epsilon)
                return Double.MaxValue;
            double crossZ = (x0 * y1 - x1 * y0) / (mag0 * mag1);
            // Z-component of the cross product of the normalized vectors
            double dot = (x0 * x1 + y0 * y1) / (mag0 * mag1); // dot product of the normalized vectors
            return Math.Atan2(crossZ, dot); // Angle between the two vectors
        }

        #endregion Contour Following
    }
}