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
            int sf = 256;
            int w = 3 * sf + 4;
            int h = 2 * sf + 4;
            var finalTransform = DenseMatrix.OfArray(new double[,] { { sf, 0, 4 * sf / 3 + 2 }, { 0, sf, sf / 3 + 2 }, { 0, 0, 1 } });
            using (DirectBitmap image = new DirectBitmap(w, h),
                image2 = new DirectBitmap(w, h))
            {
                int white = unchecked((int)0xffffffff); // white
                int black = unchecked((int)0xff000000); // black
                int red = unchecked((int)0xffff0000); // red
                int blue = unchecked((int)0xff0000ff); // blue
                //int green = unchecked((int)0xff00ff00); // green

                // Initialize the image with white
                for (int y = 0; y < h; ++y)
                    for (int x = 0; x < w; ++x)
                        image.Bits[w * y + x] = white;

                fractal.RenderFractal(image, 18, Fractal.renderLine, finalTransform, black);
                ImProc.BWBoundary(image, image2);
                fractal.ReferenceImage = image2;

                //fractal.RenderFractal(image, 17, Fractal.renderLine, finalTransform * fractal.IFS[0], green);
                fractal.RenderFractal(image, 0, Fractal.renderPrimarySpirals, finalTransform, red);
                for (int i = 0; i < 10; ++i)
                    fractal.RenderFractal(image, i, fractal.renderSecondarySpirals, finalTransform, blue);

                IO.SaveImage(image, @"C:\Users\info\Desktop\fractal.bmp");
                //IO.SaveImage(image2, @"C:\Users\info\Desktop\fractalborder.bmp");
            }
        }
    }
}
