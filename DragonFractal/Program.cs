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
            var finalTransform = DenseMatrix.OfArray(new double[,] { { 100, 0, 250 }, { 0, 100, 150 }, { 0, 0, 1 } });
            using (DirectBitmap image = new DirectBitmap(500, 500))
            {
                int white = unchecked((int)0xffffffff); // white
                for (int y = 0; y < 500; ++y)
                    for (int x = 0; x < 500; ++x)
                        image.Bits[500 * y + x] = white;
                fractal.RenderFractal(image, 16, Fractal.renderLine, finalTransform);

                int red = unchecked((int)0xffff0000); // red
                // Primary spiral
                ImProc.DrawSpiral(150, 250, Math.PI / 2, 100 * 1 / 3.0, 40 * Math.PI, Math.PI / 2, Math.PI / 512, 1.0 / 16.0, red, image);
                ImProc.DrawSpiral(350, 250, -Math.PI / 2, 100 * 2 / 3.0, 40 * Math.PI, Math.PI / 2, Math.PI / 512, 1.0 / 16.0, red, image);

                // Secondary spirals
                int blue = unchecked((int)0xff0000ff); // blue
                ImProc.DrawSpiral(250, 150, Math.PI / 2, 100 / 3.0, 40 * Math.PI, Math.PI / 4, Math.PI / 512, 1.0 / 16.0, blue, image);
                ImProc.DrawSpiral(250, 150, 3 * Math.PI / 4, 100 / 3.0, 40 * Math.PI, Math.PI / 4, Math.PI / 512, 1.0 / 16.0, blue, image);
                ImProc.DrawSpiral(250, 150, Math.PI, 100 / 3.0, 40 * Math.PI, 3 * Math.PI / 8, Math.PI / 512, 1.0 / 16.0, blue, image);

                IO.SaveImage(image, @"C:\Users\info\Desktop\fractal.bmp");
            }
        }
    }
}
