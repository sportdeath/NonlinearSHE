#include <iostream>
#include <vector>
#include <algorithm>

#include <CImg.h>

#include "imageFunctions.hpp"

void ImageFunctions::RGBToLong(long & output, 
               const unsigned char & r, 
               const unsigned char & b, 
               const unsigned char& g,
               const cimg_library::CImg<unsigned char> & colorMap) {

  // This computes the minimum Euclidean
  // distance from an RGB value to values
  // in the colorMap
  //
  // The row on the colorMap whose distance
  // from the input RGB values is returned

  double minDistance = DBL_MAX;
  output = 0;
  for (long row = 0; row < colorMap.height(); row++) {
    // Euclidean distance
    double distance = sqrt(
        pow((r - colorMap(0,row,0,0)), 2.0) +
        pow((g - colorMap(0,row,0,1)), 2.0) +
        pow((b - colorMap(0,row,0,2)), 2.0));

    if (distance < minDistance) {
      minDistance = distance;
      output = row;
    }
  }

}


void ImageFunctions::longToRGB(const long & input,
               unsigned char & r, 
               unsigned char & b, 
               unsigned char & g,
               const cimg_library::CImg<unsigned char> & colorMap) {
  r = colorMap(0,input,0,0);
  g = colorMap(0,input,0,1);
  b = colorMap(0,input,0,2);
}


void ImageFunctions::imageToLongs(
    std::vector<long> & output,
    const cimg_library::CImg<unsigned char> & image,
    const cimg_library::CImg<unsigned char> & colorMap
    ) {

  output.resize(image.height() * image.width());

  for (long r = 0; r < image.height(); r++) {
    for (long c = 0; c < image.width(); c++) {

      RGBToLong(output[c + r * image.width()], 
                image(c,r,0,0), 
                image(c,r,0,1),
                image(c,r,0,2),
                colorMap);
    }
  }
}


void ImageFunctions::longsToImage(
    const std::vector<long> & input,
    cimg_library::CImg<unsigned char> & image,
    const cimg_library::CImg<unsigned char> & colorMap
    ) {
  for (long row = 0; row < image.height(); row++) {
    for (long col = 0; col < image.width(); col++) {
      longToRGB(input[col + row * image.width()], 
                image(col,row,0,0), 
                image(col,row,0,1),
                image(col,row,0,2),
                colorMap);
    }
  }
}


/**
 * The following code is blatantly plagiarized
 * from Nan C. Schaller professor of Computer
 * Science at The Rochester Institute of Technology.
 *
 * https://www.cs.rit.edu/~ncs/color/t_convert.html
 */
void ImageFunctions::RGBtoHSV(const unsigned char rInt,
              const unsigned char gInt,
              const unsigned char bInt,
              double * h,
              double * s,
              double * v) {

  // RGB are now in [0,1)
  double r, g, b;
  r = rInt/256.;
  b = bInt/256.;
  g = gInt/256.;

  double min, max, delta;
  min = std::min(std::min(r,g),b);
  max = std::max(std::max(r,g),b);
  *v = max;
  
  delta = max - min;
  if ( max != 0) {
    *s = delta/max;
  } else {
    *s = 0;
    *h = -1;
    return;
  }

  if ( r == max) {
    *h = (g - b) /delta;
  } else if (g == max) {
    *h = 2 + (b - r)/delta;
  } else {
    *h = 4 + (r - g)/delta;
  }
  *h *= 60;

  while (*h < 0) {
    *h += 360;
  }
}



/**
 * The following code is blatantly plagiarized
 * from Nan C. Schaller professor of Computer
 * Science at The Rochester Institute of Technology.
 *
 * https://www.cs.rit.edu/~ncs/color/t_convert.html
 */
void ImageFunctions::HSVtoRGB(unsigned char * rInt,
              unsigned char * gInt,
              unsigned char * bInt,
              const double h,
              const double s,
              const double v) {

  int i;
  double f, p, q, t;

  if ( s== 0) {
    *rInt = *gInt = *bInt = v * 256.;
    return;
  }

  double fracH = h/60.;
  i = floor(fracH);
  f = fracH-i;
  p = v * (1 - s);
  q = v * (1 - s * f);
  t = v * (1 - s * (1 - f));

  double r,g,b;
  switch (i) {
    case 0:
      r = v;
      g = t;
      b = p;
      break;
    case 1:
      r = q;
      g = v;
      b = p;
      break;
    case 2:
      r = p;
      g = v;
      b = t;
      break;
    case 3:
      r = p;
      g = q;
      b = v;
      break;
    case 4:
      r = t;
      g = p;
      b = v;
      break;
    default:
      r = v;
      g = p;
      b = q;
      break;
  }
  
  *rInt = r * 256.;
  *bInt = b * 256.;
  *gInt = g * 256.;
}
