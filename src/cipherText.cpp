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


void YASHE_CT::polyRecursion(YASHE_CT & output,
                             const NTL::ZZ_pX& poly,
                             const std::vector<YASHE_CT> & powers,
                             const std::vector<YASHE_CT> & powersOfPowers) {

    std::cout << "starting... " << std::endl;

    long n = deg(poly);
    long k = powers.size() - 1;

    std::cout << "got parameters" << std::endl;

    if (n <= k) {
        // Just multiply coefficients by powers
        YASHE_CT product;
        output = YASHE_CT(0, powers[1].y);
        for (long i = 0; i <= deg(poly); i++) {
            mul(product, rem(rep(poly[i]), powers[1].y -> getPModulus()) , powers[i]);
            add(output, output, product);
        }
    } else {
        long l = std::ceil(std::log2(n/double(k)) + 1);

        while ( k * (pow(2.0, l) - 1) > n) {
          l -= 1;
        }
        while ( k * (pow(2.0, l) - 1) < n) {
          l += 1;
        }

        std::cout << "n: " << n << ", k: " << k << ", l: " << l << std::endl;

        long p = pow(2.0, l -1);

        std::cout << "p: " << p << std::endl;

        NTL::ZZ_pX q, r, c, s;
        {
          NTL::ZZ_pPush push(NTL::ZZ(powers[1].y -> getPModulus()));
          // q is quotient with x^(kp)
          // r is the remainder
          // c is quotient of r - x^(k(p - 1)) over q
          // s is the remainder
          // x^(kp)
          NTL::ZZ_pX xToKP(NTL::INIT_MONO, k * p);
          // x^(k(p -1))
          NTL::ZZ_pX xToKPMinusOne(NTL::INIT_MONO, k * (p - 1));


          std::cout << "poly: " << poly << std::endl;
          std::cout << "x^(kp): " << xToKP << std::endl;

          DivRem(q, r, poly, xToKP);
          std::cout << "quotient q: " << q << std::endl;
          std::cout << "remainder r: " << r << std::endl;
          std::cout << "x^(k(p - 1)): " << xToKPMinusOne << std::endl;
          std::cout << "r - x^k(p-1):" << r - xToKPMinusOne << std::endl;
          DivRem(c, s, r - xToKPMinusOne, q);
          std::cout << "quotient c: " << c << std::endl;
          std::cout << "remainder s: " << s << std::endl;
        }

        // q has degree k(p - 1)
        // c has degree <= k -1
        // s has degree <=k(p - 1) - 1
        YASHE_CT qOutput, cOutput, sOutput;

        polyRecursion(qOutput, q, powers, powersOfPowers);
        polyRecursion(cOutput, c, powers, powersOfPowers);
        polyRecursion(sOutput, s, powers, powersOfPowers);

        std::cout << "about to do some computation" << std::endl;
        std::cout << "addition..." << std::endl;
        std::cout << "p: " << p << std::endl;
        std::cout << powersOfPowers.size() << std::endl;
        add(output, cOutput, powersOfPowers[l - 1]);
        std::cout << "multiplication..." << std::endl;
        mul(output,  output, qOutput);
        std::cout << "addition..." << std::endl;
        add(output,  output, powersOfPowers[l - 2]);
        std::cout << "addition..." << std::endl;
        add(output,  output, sOutput);
        std::cout << "computation done" << std::endl;
    }
}


/**
 * Polynomials are evaluated using method by Paterson
 * and Stockmeyer from the 1973 paper "On the number of
 * multiplications necessary to evaluate polynomials".
 * O(sqrt(d)) ciphertext multiplications must be made
 * with a depth of O(log(d)) where d is the degree.
 */
void YASHE_CT::evalPoly(YASHE_CT& output,
                        const YASHE_CT & input,
                        NTL::ZZ_pX & poly) {
  YASHE * y = input.y;

  poly.normalize();

  long n = deg(poly);
  long k = std::floor(std::sqrt(n/2.));
  long m = std::ceil(std::log2(1. + std::sqrt(2.*n))); // n = k((2^m) - 1)

  while ( k * (pow(2.0, m) - 1) > n) {
    k -= 1;
  }
  while ( k * (pow(2.0, m) - 1) < n) {
    k += 1;
  }

  std::cout << "n: " << n << " k: " << k << " m: " << m << std::endl;

  //Compute x^2 ... x^k
  std::vector<YASHE_CT> powers(k + 1);
  std::vector<long> depth(k + 1);
  powers[0] = YASHE_CT(1, y);
  powers[1] = input;
  depth[0] = 0;
  depth[1] = 0;
  //powers[1].generateMultiplier();

  for (long i = 2; i <= k; i++) {
    std::cout << "i: " << i << std::endl;
    long minimumDepth = k;
    long minimumIndex = 0;

    for (long j = 1; j <= i/2.; j++) {
      long newDepth = std::max(depth[j],depth[i - j]);

      if (newDepth < minimumDepth) {
        minimumDepth = newDepth;
        minimumIndex = j;
      } 
    }

    depth[i] = minimumDepth + 1;

    mul(powers[i], powers[minimumIndex], powers[i - minimumIndex]);
  }

  std::cout << "computed powers" << std::endl;

  // Compute x^(2k), x^(4k), x^(8k) ... x^(2^(m-1)k)
  std::vector<YASHE_CT> powersOfPowers(m);
  powersOfPowers[0] = YASHE_CT(1, y);
  powersOfPowers[1] = powers[k];
  for (long i = 2; i < m; i++) {
      mul(powersOfPowers[i], powersOfPowers[i - 1], powersOfPowers[i - 1]);
  }

  std::cout << "computed powers of powers!" << std::endl;
  
  // recurse
  polyRecursion(output, poly, powers, powersOfPowers);

  // implimenting for sparse polynomials?? do we actually need to compute every power?
  // just a couple of the k multiplications unessesary...
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


  NTL::ZZ_pX logPoly = Functions::functionToPoly(divLog, t);
  NTL::ZZ_pX expPoly = Functions::functionToPoly(divExp, t);

  YASHE_CT exponentA, exponentB;

  evalPoly(exponentA, a, logPoly);
  evalPoly(exponentB, b, logPoly);

  sub(exponentA, exponentA, exponentB);

  evalPoly(output, exponentA, expPoly);
}
