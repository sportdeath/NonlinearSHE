#include <iostream>
#include <vector>
#include <algorithm>

#include <CImg.h>

#include <YASHE/YASHE.hpp>
#include <YASHE/cipherText.hpp>
#include <YASHE/functions.hpp>

#include "imageFunctions.hpp"

/**
 * This demonstrates homomorphic operations being done
 * on images. An input image is encrypted pixelwise and
 * homomorphically. Then a homomorphic transformation 
 * is done on it. In this case we convert the RGB values 
 * of the pixels of the image to HSV (Hue, Saturation, Value).
 * Then we rotate the hue by a constant amount, morphing
 * the color of the image. The result is then decrypted
 * and displayed next to the original.
 *
 * Due to the constraints of homomorphic computation the
 * image is scaled down to be of size 100 by 100 pixels.
 * Then we perform batch homomorphic computation on a 
 * vector of pixels of size ~10,000.
 */
int main(int argc, char * argv[]) {

  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " ImageFile1 ImageFile2" << std::endl;
    return 1;
  }

  // Our image input
  cimg_library::CImg<unsigned char> img1(argv[1]);
  cimg_library::CImg<unsigned char> img2(argv[2]);

  clock_t start, end;

  start = clock();

  // Generate parameters for the YASHE protocol
  // and create environment.
  YASHE SHE = YASHE::readFromFile("resources/8BitFHE");
  long t = SHE.getPModulus();

  end = clock();
  std::cout << "Reading parameters completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  long facts = SHE.getNumFactors();

  std::cout << facts << " factors." << std::endl;

  // Resize the image so we can pack it into a single
  // cipher text
  if (img1.width() * img1.height() > facts) {
    double scalingFactor = facts/double(img1.width() * img1.height());
    scalingFactor = sqrt(scalingFactor);
    long newWidth = img1.width() * scalingFactor;
    long newHeight = img1.height() * scalingFactor;
    img1.resize(newWidth,newHeight, 1, 3, 5);
    img2.resize(newWidth,newHeight, 1, 3, 5);
  }

  // Convert the image into a list of integers
  // each integer represents a pixel defined
  // by the color mapping above. 
  std::vector<std::vector<long>> channels1(3);
  std::vector<std::vector<long>> channels2(3);
  for (long channel = 0; channel < 3; channel ++) {
    channels1[channel].resize(facts);
    channels2[channel].resize(facts);
  }
  for (long r = 0; r < img1.height(); r++) {
    for (long c = 0; c < img1.width(); c++) {
      long index = c + r * img1.width();
      for (long channel = 0; channel < 3; channel ++) {
        channels1[channel][index] = img1(c,r,0,channel);
        channels2[channel][index] = img2(c,r,0,channel);
      }
    }
  }

  // The function is converted into
  // a polynomial of degree t = 257
  std::vector<long> divideBy2 = Functions::functionToPoly(
      [] (long input) { return input/2;}, t);

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
  std::vector<YASHE_CT> encChannels1(3);
  std::vector<YASHE_CT> encChannels2(3);
  for (long channel = 0; channel < 3; channel ++) {
    encChannels1[channel] = SHE.encryptBatch(channels1[channel]);
    encChannels2[channel] = SHE.encryptBatch(channels2[channel]);
  }

  end = clock();
  std::cout << "Encryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // evaluate the polynomial
  
  for (long channel = 0; channel < 3; channel ++) {
    YASHE_CT::evalPoly(encChannels1[channel], encChannels1[channel], divideBy2);
    YASHE_CT::evalPoly(encChannels2[channel], encChannels2[channel], divideBy2);
    YASHE_CT::add(encChannels1[channel], encChannels1[channel], encChannels2[channel]);
  }

  end = clock();
  std::cout << "Polynomial evaluation completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // decrypt the message
  std::vector<std::vector<long>> decryption(3);
  for (long channel = 0; channel < 3; channel ++) {
    decryption[channel] = SHE.decryptBatch(encChannels1[channel], secretKey);
  }

  end = clock();
  std::cout << "Decryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  // turn the message back into an image
  cimg_library::CImg<unsigned char> 
    outputImg(img1.width(), img1.height(), 1, 3);

  for (long r = 0; r < img1.height(); r++) {
    for (long c = 0; c < img1.width(); c++) {
      long index = c + r * img1.width();
      for (long channel = 0; channel < 3; channel ++) {
        outputImg(c,r,0,channel) = decryption[channel][index];
      }
    }
  }

  // Display the input next to the output!
  (img1,img2,outputImg).display("result!",false);

  return 0;
}
