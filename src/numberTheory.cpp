#include <vector>

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

#include "YASHE/numberTheory.hpp"

long NumberTheory::GCD(long a, long b) {
  return b == 0 ? a : GCD(b, a % b);
}

long NumberTheory::eulerToitient(long n) {
  long output = n;
  for (long i = 2; i * i <= n; i++) {
    // if relatively prime
    if (n % i == 0) {
      // get rid of all factors
      while (n % i == 0) {
        n /=i;
      }
      output -= output/i;
    }
  }
  if (n > 1) {
    output -= output/n;
  }
  return output;
}

/**
 * We compute the nth cyclotomic polynomial
 * by using the following recursive equation:
 *
 *                         x^n - 1
 *  Phi_n(x) = -------------------------------
 *             Product[Phi_d(x) | d divides n]
 *
 */
NTL::ZZX NumberTheory::cyclotomicPoly(long n) {

  // x^n -1
  NTL::ZZX numerator;
  NTL::SetCoeff(numerator, n);
  NTL::SetCoeff(numerator, 0, -1);

  if (n == 1) {
    return numerator;
  } else {
    NTL::ZZX denominator(1);
    for (long d = 1; 2*d <= n; d++) {
      if (n % d == 0) {
        denominator *= cyclotomicPoly(d);
      }
    }
    return numerator/denominator;
  }
}


/**
 * We compute the Chinese Remainder theorem
 * by using the equation given in the existence
 * proof:
 *  
 *  Let output = x, 
 *      inputs[i] = a_i,
 *      factors[i] = n_i, 
 *      modulus = N
 *      
 *  if x = a_i mod n_i for all i
 *  and N = Product[n_i, for all i]
 *  then
 *
 *          [       N    ( ( ( N )^-1)         )]
 *   x = Sum[a_i * --- * ( ( (---)   ) mod n_i )]
 *          [      n_i   ( ( (n_i)   )         )]
 */
void NumberTheory::CRT(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const std::vector<NTL::ZZ_pX>& factors) {

  long numFactors = factors.size();
  output = NTL::ZZ_pX::zero();

  NTL::ZZ_pXModulus modulusM(modulus);
  //NTL::ZZ_pXModulus factorMod;

  NTL::ZZ_pX fInv, fInvInv;

  for (long i = 0; i < numFactors; i++) {
    //build(factorMod, factors[i]);
    div(fInv, modulusM, factors[i]);
    rem(fInvInv, fInv, factors[i]);
    InvMod(fInvInv, fInvInv, factors[i]);

    // These values here could be computed one time and stored
    //         vvvvvvvvvvvvvvvvvvvvvvvv
    output += MulMod(fInv, fInvInv, modulusM) * inputs[i];
    // Then you could just multiply by input and sum
  }

}

void NumberTheory::CRTwithElements(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const std::vector<NTL::ZZ_pX>& crtElements) {

  output = NTL::ZZ_pX::zero();

  for (long i = 0; i < crtElements.size(); i++) {
    output += crtElements[i] * inputs[i];
  }

}
