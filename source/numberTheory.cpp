#include <vector>

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

#include "NumberTheory.hpp"

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

void NumberTheory::CRT(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const NTL::vec_ZZ_pX& factors) {

  long numFactors = factors.length();
  output = NTL::ZZ_pX::zero();

  NTL::ZZ_pXModulus modulusM(modulus);


  NTL::ZZ_pX fInv, fInvInv;

  for (long i = 0; i < numFactors; i++) {
    fInv = modulus/factors[i];

    rem(fInvInv, fInv, factors[i]);
    InvMod(fInvInv, fInvInv, factors[i]);

    output += MulMod(fInv, fInvInv, modulusM) * inputs[i];
  }

}
