#include <vector>
#include <algorithm>

#include <NTL/ZZ_pX.h>

#include "YASHE/YASHE.hpp"
#include "YASHE/cipherText.hpp"
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
  if (a.isMultiplier and b.isMultiplier) {
    a.y -> roundMultiply(output.poly, a.multiplier, b.multiplier);
  } else if (a.isMultiplier) {
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

  long t = poly.size() - 1;

  // sqrt(t + 1)
  long sqrtDegree = std::ceil(sqrt(t + 1));

  // A vector of all of the powers of the input
  // from 1 to sqrtDegree
  std::vector<YASHE_CT> powers(sqrtDegree + 1);

  powers[0] = YASHE_CT(1, y);
  powers[1] = input;
  powers[0].generateMultiplier();
  powers[1].generateMultiplier();

  for (long i = 2; i <= sqrtDegree; i++) {
    long leftIndex = std::floor(i/2.);
    long rightIndex = std::ceil(i/2.);
    mul(powers[i], powers[leftIndex], powers[rightIndex]);
  }

  // A vector of x^sqrtDegree, x^2sqrtDegree ... x^degree
  long maximumChunk = std::floor(t/double(sqrtDegree));

  std::vector<YASHE_CT> powersOfPowers(maximumChunk + 1);
  powersOfPowers[0] = YASHE_CT(1, y);
  powersOfPowers[1] = powers[sqrtDegree];
  powersOfPowers[0].generateMultiplier();
  powersOfPowers[1].generateMultiplier();

  for (long i = 2; i <= maximumChunk; i++) {
    long leftIndex = std::floor(i/2.);
    long rightIndex = std::ceil(i/2.);
    mul(powersOfPowers[i], powersOfPowers[leftIndex], powersOfPowers[rightIndex]);
  }

  // The max chunk
  output = YASHE_CT(0, y);
  for (long chunk = 0; chunk <= maximumChunk; chunk++) {

    YASHE_CT subTotal(poly[chunk * sqrtDegree], y);

    for (long i = 1; i < std::min(sqrtDegree, t - chunk*sqrtDegree + 1); i++) {
      long coefficient = poly[chunk * sqrtDegree + i];
      if (coefficient != 0) {
        YASHE_CT product;
        mul(product, poly[chunk * sqrtDegree + i], powers[i]);
        add(subTotal, subTotal, product);
      }
    }

    if (chunk != 0) {
      mul(subTotal, subTotal, powersOfPowers[chunk]);
    }
    add(output, output, subTotal);
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
  long t = a.y -> getPModulus();

  std::function<long(long)> divLog = [t](long input) {

    if (input == 0) {
      return long(0);
    } else {
      return long(std::round(t/2.*log(input)/log(t)));
    }
  };

  std::function<long(long)> divExp = [t](long input) {

    if (input > t/2.) {
      return long(0);
    } else {
      return long(std::floor(pow(t, 2.* input/ double(t))));
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

void YASHE_CT::geq(YASHE_CT & output, YASHE_CT & x, YASHE_CT & y) {
  long t = x.y -> getPModulus();

  std::function<long(long)> geqFunc = [t](long input) {
    return input >= (t - 1)/2.;
  };

  std::function<long(long)> leqFunc = [t](long input) {
    return input <= (t - 1)/2.;
  };

  std::vector<long> geqPoly = Functions::functionToPoly(geqFunc, t);
  std::vector<long> leqPoly = Functions::functionToPoly(leqFunc, t);

  YASHE_CT a, b, c, ab;

  evalPoly(a, x, geqPoly);
  evalPoly(b, y, leqPoly);

  sub(c, x, y);

  evalPoly(c, c, leqPoly);

  mul(ab, a, b);
  mul(output, -2, ab); // -2ab
  add(output, a, output); //a -2ab
  add(output, b, output); //a + b -2ab
  mul(output, c, output); // c(a+b-2ab)
  add(output, ab, output); // ab + c(a+b-2ab)
}
