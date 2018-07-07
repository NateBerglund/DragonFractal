using System;
using System.Drawing;
using System.IO;
using System.Linq;

namespace DragonFractal
{
    /// <summary>
    /// Class for reading/writing images
    /// </summary>
    static class IO
    {
        #region Image Saving

        /// <summary>
        /// Save an image
        /// </summary>
        /// <param name="image">Image to save, stored as a Bitmap</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveImage(Bitmap image, string filename)
        {
            SaveImageHelper(image, filename);
        }

        /// <summary>
        /// Save an image
        /// </summary>
        /// <param name="image">Image to save, stored as a DirectBitmap</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveImage(DirectBitmap image, string filename)
        {
            SaveImageHelper(image, filename);
        }

        /// <summary>
        /// Save an image.
        /// </summary>
        /// <param name="image">Image to save, stored as a 2-d array of bytes</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveImage(byte[,] image, string filename)
        {
            int w = image.GetLength(1);
            int h = image.GetLength(0);
            using (DirectBitmap fastImage = new DirectBitmap(w, h))
            {
                for (int y = 0; y < h; ++y)
                {
                    for (int x = 0; x < w; ++x)
                    {
                        int val = image[y, x];
                        fastImage.Bits[y * w + x] = (255 << 24) | (val << 16) | (val << 8) | val;
                    }
                }
                SaveImageHelper(fastImage, filename);
            }
        }

        /// <summary>
        /// Overlays a greyscale image with a specified weight on top of an existing color image. The greater the weight,
        /// the more the output image will resemble "overlay", and the lesser the weight, the more the output image will
        /// resemble "image".
        /// </summary>
        /// <param name="image">Original image</param>
        /// <param name="overlay">Image to overlay</param>
        /// <param name="alpha">Overlay weight between 0 and 1</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveOverlayImage(DirectBitmap image, byte[,] overlay, double alpha, string filename)
        {
            int w = image.Width;
            int h = image.Height;
            using (DirectBitmap image2 = new DirectBitmap(w, h))
            {
                for (int y = 0; y < h; ++y)
                {
                    for (int x = 0; x < w; ++x)
                    {
                        int colorValue = image.Bits[w * y + x];
                        int R = (colorValue & 0x00ff0000) >> 16;
                        int G = (colorValue & 0x0000ff00) >> 8;
                        int B = (colorValue & 0x000000ff);
                        byte V = overlay[y, x];
                        R = (byte)Math.Round((1.0 - alpha) * R + alpha * V);
                        G = (byte)Math.Round((1.0 - alpha) * G + alpha * V);
                        B = (byte)Math.Round((1.0 - alpha) * B + alpha * V);
                        image2.Bits[y * w + x] = (255 << 24) | (R << 16) | (G << 8) | B;
                    }
                }
                SaveImageHelper(image2, filename);
            }
        }

        /// <summary>
        /// Overlays an image with a specified weight on top of an existing image. The greater the weight,
        /// the more the output image will resemble "overlay", and the lesser the weight, the more the output image will
        /// resemble "image".
        /// </summary>
        /// <param name="image">Original image</param>
        /// <param name="overlay">Image to overlay</param>
        /// <param name="alpha">Overlay weight between 0 and 1</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveOverlayImage(DirectBitmap image, DirectBitmap overlay, double alpha, string filename)
        {
            int w = image.Width;
            int h = image.Height;
            using (DirectBitmap image2 = new DirectBitmap(w, h))
            {
                for (int y = 0; y < h; ++y)
                {
                    for (int x = 0; x < w; ++x)
                    {
                        int colorValue = image.Bits[w * y + x];
                        int R = (colorValue & 0x00ff0000) >> 16;
                        int G = (colorValue & 0x0000ff00) >> 8;
                        int B = (colorValue & 0x000000ff);
                        int overlayValue = overlay.Bits[w * y + x];
                        int R2 = (overlayValue & 0x00ff0000) >> 16;
                        int G2 = (overlayValue & 0x0000ff00) >> 8;
                        int B2 = (overlayValue & 0x000000ff);
                        R = (byte)Math.Round((1.0 - alpha) * R + alpha * R2);
                        G = (byte)Math.Round((1.0 - alpha) * G + alpha * G2);
                        B = (byte)Math.Round((1.0 - alpha) * B + alpha * B2);
                        image2.Bits[y * w + x] = (255 << 24) | (R << 16) | (G << 8) | B;
                    }
                }
                SaveImageHelper(image2, filename);
            }
        }

        /// <summary>
        /// Saves a comparison image show what pixels changed between two images. Pixels that were
        /// the same color in both images are shown in their original color, and those
        /// that differ are shown in the color passed in as "differenceColor".
        /// </summary>
        /// <param name="image1">First image</param>
        /// <param name="image2">Second image</param>
        /// <param name="differenceColor">Color to give to the pixels that are different (ARGB value)</param>
        /// <param name="filename">Path to save the image to</param>
        public static void SaveBWCompareImage(byte[,] image1, byte[,] image2, int differenceColor, string filename)
        {
            int w = image1.GetLength(1);
            int h = image1.GetLength(0);
            if (image2.GetLength(1) != w || image2.GetLength(0) != h)
                throw new ArgumentException("SaveBWCompareImage: input image sizes do not match!", "image2");
            using (DirectBitmap image = new DirectBitmap(w, h))
            {
                for (int y = 0; y < h; ++y)
                {
                    for (int x = 0; x < w; ++x)
                    {
                        if (image1[y, x] == image2[y, x])
                            image.Bits[w * y + x] = (255 << 24) | (image1[y, x] << 16) | (image1[y, x] << 8) | image1[y, x];
                        else
                            image.Bits[w * y + x] = differenceColor;
                    }
                }
                SaveImageHelper(image, filename);
            }
        }

        #endregion Image Saving

        #region Image Loading

        /// <summary>
        /// Load an image
        /// </summary>
        /// <param name="filename">Path to the image</param>
        /// <returns></returns>
        public static DirectBitmap LoadImage(string filename)
        {
            using (Bitmap bmp = new Bitmap(filename))
            {
                // Convert image into format that allows fast pixel access
                int w = bmp.Width;
                int h = bmp.Height;
                DirectBitmap image = new DirectBitmap(w, h);
                using (var g = Graphics.FromImage(image.Bitmap))
                {
                    g.DrawImage(bmp, 0, 0, w, h);
                }
                return image;
            }
        }

        #endregion Image Loading

        #region Private Helper Functions

        /// <summary>
        /// This is needed because image.Save(filename) sometimes causes an exception with the message
        /// "A generic error occurred in GDI+" to occur. No idea why, but this is the fix I found online for it.
        /// </summary>
        /// <param name="image">Image to save, stored as a DirectBitmap</param>
        /// <param name="filename">Path to save the image to</param>
        private static void SaveImageHelper(Bitmap image, string filename)
        {
            // To get around the GDI+ error, we save to a memory stream first, then from the stream to the file.
            using (MemoryStream memory = new MemoryStream())
            {
                using (FileStream fs = new FileStream(filename, FileMode.Create, FileAccess.ReadWrite))
                {
                    image.Save(memory, System.Drawing.Imaging.ImageFormat.Bmp);
                    byte[] bytes = memory.ToArray();
                    fs.Write(bytes, 0, bytes.Length);
                }
            }
        }

        /// <summary>
        /// This is needed because DirectBitmap stores its pixels internally as 32-bit,
        /// even though we don't actually use transparency. Some applications have
        /// trouble with 32-bit images, so this helps ensure all debug output images are 24-bit.
        /// </summary>
        /// <param name="image">Image to save, stored as a Bitmap</param>
        /// <param name="filename">Path to save the image to</param>
        private static void SaveImageHelper(DirectBitmap image, string filename)
        {
            // Since the DirectBitmap type uses 32-bit images, but we want to save a 24-bit image, we must convert the bitmap here
            using (Bitmap bm24 = new Bitmap(image.Bitmap.Width, image.Bitmap.Height, System.Drawing.Imaging.PixelFormat.Format24bppRgb))
            {
                using (Graphics g = Graphics.FromImage(bm24))
                {
                    g.DrawImage(image.Bitmap, 0, 0);
                }
                SaveImageHelper(bm24, filename);
            }
        }

        #endregion Private Helper Functions
    }
}
