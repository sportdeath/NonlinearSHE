#include <vector>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include "functions.hpp"

std::vector<long> Functions::functionToPoly(std::function<long(long)> f, long modulus) {
  // Temporarily switch the modulus to the new modulus
  NTL::ZZ bigModulus(modulus);
  NTL::ZZ_pPush push(bigModulus);

  // Make the polynomial modulus x^modulus - x
  // By Fermat's Little Theorem x^modulus - x = 0
  // (modulus must be prime)
  NTL::ZZ_pX polyModulusX;
  NTL::SetCoeff(polyModulusX, modulus);
  NTL::SetCoeff(polyModulusX, 1, -1);
  NTL::ZZ_pXModulus polyModulus(polyModulusX);

  // Zero polynomial
  NTL::ZZ_pX out;
  out.SetLength(modulus);

  // The delta function
  // 1 if x = i, 0 otherwise
  // delta(x - i) = 1 - (x - i)^(p - 1)
  NTL::ZZ_pX deltaX;

  // out += f(i)*delta(x - i)
  for (long i = 0; i < modulus; i++) {
    deltaX = PowerXPlusAMod(NTL::ZZ_p(-i), modulus - 1, polyModulus);
    deltaX = 1 - deltaX;
    out += f(i) * deltaX;
  }

  std::vector<long> vectorOut(modulus);

  for (long i = 0; i < modulus; i++) {
    vectorOut[i] = rem(rep(out[i]), modulus);
  }

  return vectorOut;
}


//std::vector<long> Functions::functionToPoly2D(
    //std::function<long(long, long)> f, long modulus) {
  //// Temporarily switch the modulus to the new modulus
  //NTL::ZZ bigModulus(modulus);
  //NTL::ZZ_pPush push(bigModulus);

  //// Make the polynomial modulus x^modulus - x
  //// By Fermat's Little Theorem x^modulus - x = 0
  //// (modulus must be prime)
  //NTL::ZZ_pX polyModulusX;
  //NTL::SetCoeff(polyModulusX, modulus);
  //NTL::SetCoeff(polyModulusX, 1, -1);
  //NTL::ZZ_pXModulus polyModulus(polyModulusX);

  //// Zero polynomial
  //NTL::ZZ_pX out;
  //out.SetLength(modulus);

  //// The delta function
  //// 1 if x = i, 0 otherwise
  //// delta(x - i) = 1 - (x - i)^(p - 1)
  //NTL::ZZ_pX deltaX;

  //// out += f(i)*delta(x - i)
  //for (long i = 0; i < modulus; i++) {
    //deltaX = PowerXPlusAMod(NTL::ZZ_p(-i), modulus - 1, polyModulus);
    //deltaX = 1 - deltaX;
    //out += f(i) * deltaX;
  //}

  //std::vector<long> vectorOut(modulus);

  //for (long i = 0; i < modulus; i++) {
    //vectorOut[i] = rem(rep(out[i]), modulus);
  //}

  //return vectorOut;
//}


std::function<long(long)> Functions::divideByConstant(long denominator) {
  return [denominator](long x) {return x/denominator;};
}

std::function<long(long)> Functions::gt(long rightSide) {
  return [rightSide](long x) {return x > rightSide;};
}

std::function<long(long)> Functions::lt(long rightSide) {
  return [rightSide](long x) {return x < rightSide;};
}

std::function<long(long)> Functions::eq(long rightSide) {
  return [rightSide](long x) {return x == rightSide;};
}

std::function<long(long)> Functions::leq(long rightSide) {
  return [rightSide](long x) {return x <= rightSide;};
}

std::function<long(long)> Functions::geq(long rightSide) {
  return [rightSide](long x) {return x >= rightSide;};
}
