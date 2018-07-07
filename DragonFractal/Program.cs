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
                fractal.RenderFractal(image, 20, Fractal.renderLine, finalTransform);
                IO.SaveImage(image, @"C:\Users\info\Desktop\fractal.bmp");
            }
        }
    }
}
