#include <vector>
#include <algorithm>

#include <NTL/ZZ_pX.h>

#include "yashe.hpp"
#include "cipherText.hpp"
#include "functions.hpp"

YASHE_CT::YASHE_CT() {
  poly = NTL::ZZ_pX(0);
  isMultiplier = false;
}


YASHE_CT::YASHE_CT(long value, YASHE * y_) {
  y = y_;
  poly = NTL::ZZ_pX(value * y -> getModulusRatio());
  isMultiplier = false;
}


YASHE_CT::YASHE_CT(NTL::ZZ_pX poly_, YASHE * y_) {
  poly = poly_;
  y = y_;
  isMultiplier = false;
}


NTL::ZZ_pX YASHE_CT::getPoly() {
  return poly;
}


/**
 * Addition of plaintexts is addition of ciphertexts
 */
void YASHE_CT::add(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b) {
  NTL::add(output.poly, a.poly, b.poly);
  output.y = a.y;
  output.isMultiplier = false;
}
void YASHE_CT::sub(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b) {
  NTL::sub(output.poly, a.poly, b.poly);
  output.y = a.y;
  output.isMultiplier = false;
}

void YASHE_CT::generateMultiplier() {
  NTL::ZZ_pPush push(y -> getBigModulus());
  multiplier = NTL::ZZ_pXMultiplier(poly, y -> getBigCycloMod());
  this -> isMultiplier = true;
}



/**
 * In order to multiply ciphertexts, we perform multiplication
 * with rounding and then a key switching procedure.
 */
void YASHE_CT::mul(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b) {
  if (a.isMultiplier) {
    a.y -> roundMultiply(output.poly, b.poly, a.multiplier);
  } else if (b.isMultiplier) {
    a.y -> roundMultiply(output.poly, a.poly, b.multiplier);
  } else {
    a.y -> roundMultiply(output.poly, a.poly, b.poly);
  }
  a.y -> keySwitch(output.poly, output.poly);
  output.y = a.y;
  output.isMultiplier = false;
}


/**
 * In order to multiply by a constant we simply need to
 * multiply by that constant
 */
void YASHE_CT::mul(YASHE_CT& output, const YASHE_CT& a, const long& b) {
  output.poly = a.poly * b;
  output.y = a.y;
  output.isMultiplier = false;
}
void YASHE_CT::mul(YASHE_CT& output, const long& a, const YASHE_CT& b) {
  mul(output, b, a);
}


/**
 * In order to add a plaintext to a constant we must multiply
 * it by the modulus ratio q/t
 */
void YASHE_CT::add(YASHE_CT& output, const YASHE_CT& a, const long& b) {
  NTL::add(output.poly, a.poly, NTL::ZZ_pX(b* a.y -> getModulusRatio()));
  output.y = a.y;
  output.isMultiplier = false;
}
void YASHE_CT::add(YASHE_CT& output, const long& a, const YASHE_CT& b) {
  add(output, b, a);
}


/**
 * Evaluated using method by Paterson and Stockmeyer.
 * Greatly reduced number of multiplications nessesary.
 */
void YASHE_CT::evalPoly(YASHE_CT& output,
                        YASHE_CT& input,
                        const std::vector<long>& poly) {

  YASHE * y = input.y;

  long sqrtDegree = sqrt(poly.size() - 1);

  while ( (sqrtDegree + 1) * sqrtDegree < poly.size()) {
    sqrtDegree += 1;
  }

  // A vector of all of the powers of the input
  // from 1 to sqrtDegree
  std::vector<YASHE_CT> powers(sqrtDegree);
  powers[0] = input;
  powers[0].generateMultiplier();
  for (long i = 1; i < sqrtDegree; i++) {
    if ( i + 1 % 2 == 0) {
      mul(powers[i], powers[(i + 1)/2 - 1], powers[(i + 1)/2 - 1]);
    } else {
      mul(powers[i], powers[i - 1], powers[0]);
    }
  }

  // A vector of x^sqrtDegree, x^2sqrtDegree ... x^degree
  std::vector<YASHE_CT> powersOfPowers(sqrtDegree);
  powersOfPowers[0] = powers[sqrtDegree - 1];
  powersOfPowers[0].generateMultiplier();
  for (long i = 1; i < sqrtDegree; i++) {
    if ( i + 1 % 2 == 0) {
      mul(powersOfPowers[i], 
          powersOfPowers[(i + 1)/2 - 1], 
          powersOfPowers[(i + 1)/2 - 1]);
    } else {
      mul(powersOfPowers[i], powersOfPowers[i - 1], powersOfPowers[0]);
    }
  }

  // Chunk 0 
  output = YASHE_CT(poly[0], y);
  YASHE_CT product;
  for (long i = 0; i < std::min(sqrtDegree - 1, long(poly.size()) - 1); i++) {
    mul(product, poly[i + 1], powers[i]);
    add(output, output, product);
  }

  long maxChunk = sqrtDegree;
  while (poly.size() - 1 < sqrtDegree * maxChunk) {
    maxChunk -= 1;
  }

  for (long chunk = 0; chunk < maxChunk; chunk++) {
    YASHE_CT subTotal(poly[(chunk + 1) * sqrtDegree], y);
    for (long i = 0; 
         i < std::min(sqrtDegree - 1, long(poly.size()) - 1 - (chunk + 1) * sqrtDegree);
         i++) {
      mul(product, poly[i + 1 + (chunk + 1) * sqrtDegree], powers[i]);
      add(subTotal, subTotal, product);
    }
    mul(product, subTotal, powersOfPowers[chunk]);
    add(output, output, product);
  }
}

void YASHE_CT::div(YASHE_CT& output, YASHE_CT& a, YASHE_CT& b) {
  // computes the log base 2 of the input.
  // Only need log_2log_2 p integer bits
  // log_2 p - log_2log_2 p - 1 floating bits
  // 1 bit which will become high if the input is less than the output

  long t = a.y -> getPModulus();
  double base = 2.02;

  std::function<long(long)> divLog = [t, base](long input) {

    if (input == 0) {
      return long(0);
    }

    double logInput = log(input)/log(base);

    double shiftAmmount = log(t)/log(base) - log(log(t)/log(base))/log(base) - 1;

    long output = round(logInput * pow(base, shiftAmmount));

    return output;
  };


  std::function<long(long)> divExp = [t, base](long input) {
    double maxInt = log(t)/log(base);
    double shiftAmmount = log(t)/log(base) - log(log(t)/log(base))/log(base) - 1;
    double actualValue = input / pow(base, shiftAmmount);

    if (actualValue > maxInt) {
      return long(0);
    } else {
      return long(pow(base, actualValue));
    }
  };


  std::vector<long> logPoly = Functions::functionToPoly(divLog, t);
  std::vector<long> expPoly = Functions::functionToPoly(divExp, t);

  YASHE_CT exponentA, exponentB;

  evalPoly(exponentA, a, logPoly);
  evalPoly(exponentB, b, logPoly);

  sub(exponentA, exponentA, exponentB);

  evalPoly(output, exponentA, expPoly);
}
