using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Runtime.InteropServices;

namespace DragonFractal
{
    /// <summary>
    /// Class that represents a bitmap and allows fast direct access to the pixels.
    /// For more details, see http://stackoverflow.com/questions/24701703/c-sharp-faster-alternatives-to-setpixel-and-getpixel-for-bitmaps-for-windows-f.
    /// This class allows us to implement image transforms (such as the "unrolling" of the circular print) in an efficient manner.
    /// Note that it is currently still necessary to read the data into some other object such as a Bitmap and then copy the data
    /// into an instance of this object (TODO: implement a constructor that reads the image directly from file).
    /// </summary>
    public class DirectBitmap : IDisposable
    {
        /// <summary>
        /// Allows access to standard Bitmap members.
        /// </summary>
        public Bitmap Bitmap { get; private set; }

        /// <summary>
        /// The bitmap data, formatted as PARGB (premultiplied alpha + RGB)
        /// </summary>
        public Int32[] Bits { get; private set; }

        /// <summary>
        /// Tracks whether the object has been disposed or not
        /// </summary>
        public bool Disposed { get; private set; }

        /// <summary>
        /// Height of the bitmap in pixels
        /// </summary>
        public int Height { get; private set; }

        /// <summary>
        /// Width of the bitmap in pixels
        /// </summary>
        public int Width { get; private set; }

        /// <summary>
        /// Handle to the data in memory
        /// </summary>
        protected GCHandle BitsHandle { get; private set; }

        /// <summary>
        /// Constructs a DirectBitmap of the desired size. Caller must Dispose this object when finished with it.
        /// </summary>
        /// <param name="width">Width of the bitmap in pixels</param>
        /// <param name="height">Height of the bitmap in pixels</param>
        public DirectBitmap(int width, int height)
        {
            Width = width;
            Height = height;
            Bits = new Int32[width * height];
            BitsHandle = GCHandle.Alloc(Bits, GCHandleType.Pinned);
            Bitmap = new Bitmap(width, height, width * 4, PixelFormat.Format32bppPArgb, BitsHandle.AddrOfPinnedObject());
        }

        /// <summary>
        /// Dispose this object, freeing the memory associated with it.
        /// </summary>
        public void Dispose()
        {
            if (Disposed) return;
            Disposed = true;
            Bitmap.Dispose();
            BitsHandle.Free();
        }
    }
}
