#include <vector>

#include <CImg.h>

class ImageFunctions {
  public:
    static void RGBToLong(long & output, 
                          const unsigned char & r, 
                          const unsigned char & b, 
                          const unsigned char& g,
                          const cimg_library::CImg<unsigned char> & colorMap);

    static void longToRGB(const long & input,
                          unsigned char & r, 
                          unsigned char & b, 
                          unsigned char & g,
                          const cimg_library::CImg<unsigned char> & colorMap);


    static void imageToLongs(
        std::vector<long> & output,
        const cimg_library::CImg<unsigned char> & image,
        const cimg_library::CImg<unsigned char> & colorMap);


    static void longsToImage(
        const std::vector<long> & input,
        cimg_library::CImg<unsigned char> & image,
        const cimg_library::CImg<unsigned char> & colorMap
        );


    static void RGBtoHSV(const unsigned char rInt,
                         const unsigned char gInt,
                         const unsigned char bInt,
                         double * h,
                         double * s,
                         double * v);


    static void HSVtoRGB(unsigned char * rInt,
                         unsigned char * gInt,
                         unsigned char * bInt,
                         const double h,
                         const double s,
                         const double v);
};
