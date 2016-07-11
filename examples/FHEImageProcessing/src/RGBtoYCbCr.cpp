#include <iostream>
#include <vector>
#include <algorithm>

#include <CImg.h>

#include <YASHE/YASHE.hpp>
#include <YASHE/cipherText.hpp>
#include <YASHE/functions.hpp>

#include <NTL/ZZ_pX.h>

#include "imageFunctions.hpp"

std::function<long(long)> multiplyByConstant(double c) {
  return [c] (long input) {
    return c * input;
  };
};

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

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " ImageFile" << std::endl;
    return 1;
  }

  // Our image input
  cimg_library::CImg<unsigned char> img(argv[1]);

  clock_t start, end;

  start = clock();

  // Generate parameters for the YASHE protocol
  // and create environment.
  YASHE SHE = YASHE::readFromFile("resources/8BitFHE");

  end = clock();
  std::cout << "Reading parameters completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  std::cout << SHE.getNumFactors() << " factors." << std::endl;

  // Resize the image so we can pack it into a single
  // cipher text
  if (img.width() * img.height() > SHE.getNumFactors()) {
    double scalingFactor = SHE.getNumFactors()/double(img.width() * img.height());
    scalingFactor = sqrt(scalingFactor);
    long newWidth = img.width() * scalingFactor;
    long newHeight = img.height() * scalingFactor;
    img.resize(newWidth,newHeight, 1, 3, 5);
  }

  // Convert the image into a list of integers
  // each integer represents a pixel defined
  // by the color mapping above. 
  std::vector<long> R(SHE.getNumFactors()), G(SHE.getNumFactors()), B(SHE.getNumFactors());
  for (long r = 0; r < img.height(); r++) {
    for (long c = 0; c < img.width(); c++) {
      long index = c + r * img.width();
      R[index] = img(c,r,0,0);
      G[index] = img(c,r,0,1);
      B[index] = img(c,r,0,2);
    }
  }

  // The function is converted into
  // a polynomial of degree t = 257
  NTL::ZZ_pX YR = Functions::functionToPoly(multiplyByConstant(0.299), 257);
  NTL::ZZ_pX YG = Functions::functionToPoly(multiplyByConstant(0.587), 257);
  NTL::ZZ_pX YB = Functions::functionToPoly(multiplyByConstant(0.114), 257);
  NTL::ZZ_pX CbR = Functions::functionToPoly(multiplyByConstant(-0.169), 257);
  NTL::ZZ_pX CbG = Functions::functionToPoly(multiplyByConstant(-0.331), 257);
  NTL::ZZ_pX CbB = Functions::functionToPoly(multiplyByConstant(0.500), 257);
  NTL::ZZ_pX CrR = Functions::functionToPoly(multiplyByConstant(0.500), 257);
  NTL::ZZ_pX CrG = Functions::functionToPoly(multiplyByConstant(-0.419), 257);
  NTL::ZZ_pX CrB = Functions::functionToPoly(multiplyByConstant(-0.081), 257);

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
  YASHE_CT cR = SHE.encryptBatch(R);
  YASHE_CT cG = SHE.encryptBatch(G);
  YASHE_CT cB = SHE.encryptBatch(B);

  end = clock();
  std::cout << "Encryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // evaluate the polynomial
  YASHE_CT cYR, cYG, cYB, cCbR, cCbG, cCbB, cCrR, cCrG, cCrB;
  YASHE_CT::evalPoly(cYR, cR, YR);
  YASHE_CT::evalPoly(cYG, cG, YG);
  YASHE_CT::evalPoly(cYB, cB, YB);
  YASHE_CT::evalPoly(cCbR, cR, CbR);
  YASHE_CT::evalPoly(cCbG, cG, CbG);
  YASHE_CT::evalPoly(cCbB, cB, CbB);
  YASHE_CT::evalPoly(cCrR, cR, CrR);
  YASHE_CT::evalPoly(cCrG, cG, CrG);
  YASHE_CT::evalPoly(cCrB, cB, CrB);

  end = clock();
  std::cout << "Polynomial evaluation completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  YASHE_CT cY, cCb, cCr;
  YASHE_CT::add(cY, cYR, cYG);
  YASHE_CT::add(cY, cY, cYB);

  YASHE_CT::add(cCb, cCbR, cCbG);
  YASHE_CT::add(cCb, cCb, cCbB);
  YASHE_CT::add(cCb, cCb, 128);

  YASHE_CT::add(cCr, cCrR, cCrG);
  YASHE_CT::add(cCr, cCr, cCrB);
  YASHE_CT::add(cCr, cCr, 128);

  end = clock();
  std::cout << "Addition completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;
  start = clock();

  // decrypt the message
  std::vector<long> Y = SHE.decryptBatch(cY, secretKey);
  std::vector<long> Cb = SHE.decryptBatch(cCb, secretKey);
  std::vector<long> Cr = SHE.decryptBatch(cCr, secretKey);

  end = clock();
  std::cout << "Decryption completed in "
            << double(end - start)/CLOCKS_PER_SEC
            << " seconds"
            << std::endl;

  // turn the message back into an image
  cimg_library::CImg<unsigned char> 
    YImg,
    CbImg,
    CrImg,
    outputImg(img.width(), img.height(), 1, 3);
  YImg.fill(255,255,255);
  CbImg.fill(255,255,255);
  CrImg.fill(255,255,255);
  for (long r = 0; r < img.height(); r++) {
    for (long c = 0; c < img.width(); c++) {
      long index = c + r * img.width();
      outputImg(c,r,0,0) = Y[index];
      outputImg(c,r,0,1) = Cb[index];
      outputImg(c,r,0,2) = Cr[index];
    }
  }
  YImg = outputImg.get_channel(0);
  CbImg = outputImg.get_channel(1);
  CrImg = outputImg.get_channel(2);
  outputImg.YCbCrtoRGB();

  // Display the input next to the output!
  (img,YImg,CbImg,CrImg,outputImg).display("result!",false);

  return 0;
}
