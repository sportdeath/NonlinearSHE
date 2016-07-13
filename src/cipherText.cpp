#include <vector>
#include <algorithm>

#include <NTL/ZZ_pX.h>

#include "YASHE/YASHE.hpp"
#include "YASHE/ciphertext.hpp"
#include "YASHE/functions.hpp"

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


void YASHE_CT::generateMultiplier() {
  NTL::ZZ_pPush push(y -> getBigModulus());
  multiplier = NTL::ZZ_pXMultiplier(poly, y -> getBigCycloMod());
  this -> isMultiplier = true;
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


/**
 * In order to multiply ciphertexts, we perform multiplication
 * with rounding and then a key switching procedure.
 */
void YASHE_CT::mul(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b) {
  // If the multiplier exists, performing multiplication
  // is much faster.
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
 * Polynomials are evaluated using method by Paterson
 * and Stockmeyer from the 1973 paper "On the number of
 * multiplications necessary to evaluate polynomials".
 * O(sqrt(d)) ciphertext multiplications must be made
 * with a depth of O(log(d)) where d is the degree.
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
  std::vector<long> depth(sqrtDegree);
  powers[0] = input;
  depth[0] = 0;
  powers[0].generateMultiplier();
  for (long i = 1; i < sqrtDegree; i++) {
    long minimumDepth = sqrtDegree;
    long secondaryDepth = sqrtDegree;
    long minimumIndex = 0;
    for (long j = 0; j < i/2 + 1; j++) {
      long newDepth = std::max(depth[j],depth[i - j - 1]);
      long newSecDepth = std::min(depth[j],depth[i - j - 1]);
      if (newDepth < minimumDepth) {
        minimumDepth = newDepth;
        minimumIndex = j;
        secondaryDepth = newSecDepth;
      } else if ( newDepth == minimumDepth) {
        if (newSecDepth < secondaryDepth) {
          minimumIndex = j;
          secondaryDepth = newSecDepth;
        }
      }
    }
    depth[i] = minimumDepth + 1;
    mul(powers[i], powers[minimumIndex], powers[i - minimumIndex -1]);
  }

  // A vector of x^sqrtDegree, x^2sqrtDegree ... x^degree
  std::vector<YASHE_CT> powersOfPowers(sqrtDegree);
  powersOfPowers[0] = powers[sqrtDegree - 1];
  depth[0] = 0;
  powersOfPowers[0].generateMultiplier();
  for (long i = 1; i < sqrtDegree; i++) {
    long minimumDepth = sqrtDegree;
    long secondaryDepth = sqrtDegree;
    long minimumIndex = 0;
    for (long j = 0; j < i/2 + 1; j++) {
      long newDepth = std::max(depth[j],depth[i - j - 1]);
      long newSecDepth = std::min(depth[j],depth[i - j - 1]);
      if (newDepth < minimumDepth) {
        minimumDepth = newDepth;
        minimumIndex = j;
        secondaryDepth = newSecDepth;
      } else if ( newDepth == minimumDepth) {
        if (newSecDepth < secondaryDepth) {
          minimumIndex = j;
          secondaryDepth = newSecDepth;
        }
      }
    }
    depth[i] = minimumDepth + 1;
    mul(powersOfPowers[i], powersOfPowers[minimumIndex], powersOfPowers[i - minimumIndex - 1]);
  }

  // The first chunk
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

  // The other chunks
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

/**
 * Division is computed by taking note of the following equality
 *
 *         a/b = 2^(log2(a) - log2(b))
 *
 * This transforms the division problem, which is a function of
 * two variables, into a linear combination of functions of
 * single variables - significantly reducing the computation
 * time for homomorphic operations. Exponentiation and logarithms
 * are just as hard as evaluating any other nonlinear function
 * of a single variable. Because the number of floating bits
 * is limited, the reduction is approximate.
 */
void YASHE_CT::div(YASHE_CT& output, YASHE_CT& a, YASHE_CT& b) {
  // When we compute the log of an input, we will
  // need no more than log2(log2(p)) bits for integers
  // 1 bit will be used to determine the denominator is
  // greater than the numerator (preventing the exponent from
  // wrapping modulo p). The rest - log2(p) - log2(log2(p)) - 1
  // bits are for floating point operations. 
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
