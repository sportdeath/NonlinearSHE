#include <vector>

#include <CImg.h>

class ImageFunctions {
  public:
    /**
     * Outputs a long representing the RGB value
     * in the colorMap that is closest to the RGB
     * value presented in the input.
     */
    static void RGBToLong(long & output, 
                          const unsigned char & r, 
                          const unsigned char & b, 
                          const unsigned char& g,
                          const cimg_library::CImg<unsigned char> & colorMap);

    /**
     * Outputs the RGB vector present at row
     * "input" of the colorMap.
     */
    static void longToRGB(const long & input,
                          unsigned char & r, 
                          unsigned char & b, 
                          unsigned char & g,
                          const cimg_library::CImg<unsigned char> & colorMap);

    /**
     * Converts an image into a vector of
     * longs by assigning each pixel a value
     * associated with colorMap.
     *
     * The vector represents the image unrolled
     * rowwise. Therefore the entry c + r*w
     * of the vector represents the pixel at
     * row r and column c in an image of width
     * w.
     */
    static void imageToLongs(
        std::vector<long> & output,
        const cimg_library::CImg<unsigned char> & image,
        const cimg_library::CImg<unsigned char> & colorMap);


    /**
     * Converts a vector of longs representing
     * pixel values into an image. This is a pseudo
     * inverse to the imageToLongs function. A certain
     * amount of quantization is performed when
     * assigning pixels to color values on the colorMap
     * which is irreversible. 
     *
     * The width and height information are taken
     * from the image, so the image must be resized
     * beforehand. The incoming vector is taken to
     * be unrolled rowwise. Therefore the entry
     * c + r*w of the vector will become the pixel
     * at row r and column c in the resulting image
     * of width w.
     */
    static void longsToImage(
        const std::vector<long> & input,
        cimg_library::CImg<unsigned char> & image,
        const cimg_library::CImg<unsigned char> & colorMap
        );


    /**
     * Converts RGB values in the range [0,256)
     * into HSV values where H is in [0,360)
     * and S and V are in [0,1]. This is the inverse
     * of the HSVtoRGB function.
     */
    static void RGBtoHSV(const unsigned char rInt,
                         const unsigned char gInt,
                         const unsigned char bInt,
                         double * h,
                         double * s,
                         double * v);


    /**
     * Returns RGB values in the range [0,256)
     * from HSV values where H is in [0,360)
     * and S and V are in [0,1]. This is the
     * inverse of the RGBtoHSV function.
     */
    static void HSVtoRGB(unsigned char * rInt,
                         unsigned char * gInt,
                         unsigned char * bInt,
                         const double h,
                         const double s,
                         const double v);
};
