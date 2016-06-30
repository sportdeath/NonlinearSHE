#include <vector>

#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include "../source/numberTheory.hpp"

bool testGCD() {
  return NumberTheory::GCD(54, 24) == 6;
}

bool testEuler() {
  
  std::vector<long> values = {   1, 1, 2, 2, 4, 2, 6, 4, 6,
                              4,10, 4,12, 6, 8, 8,16, 6,18,
                              8,12,10,22, 8,20,12,18,12,28};

  bool areSame = true;

  for (long i = 0; i < values.size(); i++) {
    areSame &= (values[i] == NumberTheory::eulerToitient(i + 1));
  }

  return areSame;
}

bool testCyclotomic() {
  NTL::ZZ_p::init(NTL::ZZ(7));

  std::vector<long> out1 { -1,  1};            // x - 1
  std::vector<long> out2 {  1,  1};            // x + 1
  std::vector<long> out3 {  1,  1,  1};        // x^2 + x + 1
  std::vector<long> out4 {  1,  0,  1};        // x^2 + 1
  std::vector<long> out5 {  1,  1,  1,  1,  1};// x^4 + x^3 + x^2 + x^1
  std::vector<long> out6 {  1, -1,  1};        // x^2 - x + 1
  std::vector<long> out7 {  1,  1,  1,  1,  1,  1,  1};
  std::vector<long> out8 {  1,  0,  0,  0,  1};
  std::vector<long> out9 {  1,  0,  0,  1,  0,  0,  1};
  std::vector<long> ou10 {  1, -1,  1, -1,  1};
  std::vector<long> ou11 {  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1};
  std::vector<long> ou12 {  1,  0, -1,  0,  1};

  std::vector< std::vector<long> > outs = {out1,
                                           out2,
                                           out3,
                                           out4,
                                           out5,
                                           out6,
                                           out7,
                                           out8,
                                           out9,
                                           ou10,
                                           ou11,
                                           ou12};

  bool areSame = true;

  for (long i = 0; i < outs.size(); i++) {
    NTL::ZZX poly = NumberTheory::cyclotomicPoly(i + 1);
    for (long j = 0; j <= NumberTheory::eulerToitient(i + 1); j++) {
      areSame &= (poly[j] == outs[i][j]);
    }
  }

  return areSame;
}

bool testCRT() {
  long q = 257;
  NTL::ZZ_p::init(NTL::ZZ(q));

  srand(time(0));

  NTL::ZZX cyclo = NumberTheory::cyclotomicPoly(16);
  NTL::ZZ_pX cycloMod = NTL::conv<NTL::ZZ_pX>(cyclo);

  NTL::vec_ZZ_pX factors = NTL::SFBerlekamp(cycloMod);

  long numFactors = factors.length();

  std::vector<long> inputs(numFactors);
  for (long i = 0; i < numFactors; i++) {
    inputs[i] = rand() % q;
  }

  
  NTL::ZZ_pXModulus modulus(cycloMod);

  // Compute CRT
  NTL::ZZ_pX output;
  NumberTheory::CRT(output, inputs, modulus, factors);

  bool isRight = true;
  // make sure that crt mod f[i] = m[i]
  
  NTL::ZZ_pX remainder;
  for (long i = 0; i < numFactors; i++) {
    NTL::rem(remainder, output, factors[i]);
    isRight &= NTL::rem(rep(remainder[0]), q) == inputs[i];
  }

  return isRight;
}

int main() {
  if (not testGCD() ) {
    std::cout << "GCD test failed!" << std::endl;
  }
  if (not testEuler() ) {
    std::cout << "Euler test failed!" << std::endl;
  }
  if (not testCyclotomic() ) {
    std::cout << "Cyclotomic poly test failed!" << std::endl;
  }
  if (not testCRT() ) {
    std::cout << "Chinese Remainder Theorem test failed!" << std::endl;
  }
}
