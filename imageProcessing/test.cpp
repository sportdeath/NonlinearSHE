#include <iostream>
#include <vector>
#include <algorithm>

#include <CImg.h>

#include "../source/yashe.hpp"
#include "../source/cipherText.hpp"
#include "../source/functions.hpp"

void RGBToLong(long & output, 
               const unsigned char & r, 
               const unsigned char & b, 
               const unsigned char& g,
               const cimg_library::CImg<unsigned char> & colorMap) {

  double minDistance = 1000;
  output = 0;
  for (long row = 0; row < colorMap.height(); row++) {
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


void longToRGB(unsigned char & r, 
               unsigned char & b, 
               unsigned char & g,
               const long & input,
               const cimg_library::CImg<unsigned char> & colorMap) {

  r = colorMap(0,input,0,0);
  g = colorMap(0,input,0,1);
  b = colorMap(0,input,0,2);
}


void imageToLongs(
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


void longsToImage(
    cimg_library::CImg<unsigned char> & image,
    const std::vector<long> & input,
    const cimg_library::CImg<unsigned char> & colorMap
    ) {
  for (long row = 0; row < image.height(); row++) {
    for (long col = 0; col < image.width(); col++) {
      longToRGB(image(col,row,0,0), 
                image(col,row,0,1),
                image(col,row,0,2),
                input[col + row * image.width()], 
                colorMap);
    }
  }
}


void RGBtoHSV(const unsigned char rInt,
              const unsigned char gInt,
              const unsigned char bInt,
              double * h,
              double * s,
              double * v) {

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


void HSVtoRGB(unsigned char * rInt,
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



int main() {
  // Our image input
  cimg_library::CImg<unsigned char> img("../resources/test2.jpg");

  clock_t start, end;

  start = clock();

  // Generate parameters for the YASHE protocol
  // and create environment.
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(1600);
  long d = 66048; // 2^9*3*43 - 1075 irreducible factors
  long sigma = 8;
  NTL::ZZ w = NTL::power2_ZZ(300);
  YASHE SHE(t,q,d,sigma,w);

  end = clock();
  std::cout << "Factorization completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  // Resize the image so we can pack it into a single
  // cipher text
  long imgDim = sqrt(SHE.getNumFactors());
  img.resize(imgDim,imgDim, 1, 3, 5);

  // Define a color space of 256 colors
  cimg_library::CImg<unsigned char> colorMap =
    cimg_library::CImg<unsigned char>::default_LUT256();

  // Convert the image into a list of integers
  // each integer represents a pixel defined
  // by the color mapping above. 
  std::vector<long> message;
  imageToLongs(message, img, colorMap);

  // Resize the message so that it fills the entire
  // message space (the end of it will be junk)
  message.resize(SHE.getNumFactors());

  // Define a function on pixels.
  // This function takes a pixel value (0 - 255)
  // and returns another pixel value representing
  // the inversion of that pixel value.
  std::function<long(long)> invertColors = [colorMap](long input) {
    unsigned char r, g, b;
    double h, s, v;
    longToRGB(r, g, b, input, colorMap);
    RGBtoHSV(r, g, b, &h, &s, &v);

    // invert colors
    h = fmod(h + 30, 360);

    HSVtoRGB(&r, &g, &b, h, s, v);
    RGBToLong(input, r, g, b, colorMap);

    return input;
  };
  // The function is converted into
  // a polynomial of degree t = 257
  std::vector<long> poly = Functions::functionToPoly(invertColors, t);

  start = clock();

  // Generate public, secret and evaluation keys
  NTL::ZZ_pX secretKey = SHE.keyGen();

  end = clock();
  std::cout << "Key generation completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // encrypt the message
  YASHE_CT ciphertext = SHE.encryptBatch(message);

  end = clock();
  std::cout << "Encryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // evaluate the polynomial
  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);

  end = clock();
  std::cout << "Polynomial evaluation completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // decrypt the message
  std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

  end = clock();
  std::cout << "Decryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  // turn the message back into an image
  cimg_library::CImg<unsigned char> outputImg(imgDim, imgDim, 1, 3);
  longsToImage(outputImg, decryption, colorMap);

  // Display the input next to the output!
  (img, outputImg).display("result!");

  return 0;
}


