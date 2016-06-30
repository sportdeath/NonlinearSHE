#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>

#include "Yashe.h"
#include "test_yashe.h"

//NTL::ZZ_p (*f)(NTL::ZZ_p) intToModFunction(long (*f)(long)) {
//}


NTL::ZZ_pX functionToCoefficients(
    NTL::ZZ_p (*f)(NTL::ZZ_p)
    ) {

  // due to the constraints of ZZ_pX
  if (NTL::ZZ_p::modulus() > std::numeric_limits<long>::max()) {
    std::cout << "modulus too large!";
  }

  long modulus = NTL::ZZ_p::modulus() % std::numeric_limits<long>::max();

  // polyMod = x^p - x (By Fermat's Little Theorem)
  NTL::ZZ_pX x;
  SetX(x);
  NTL::ZZ_pXModulus polyModulus(power(x, modulus) - x);

  // Zero polynomial
  NTL::ZZ_pX out;

  NTL::ZZ_p index;
  NTL::ZZ_pX deltaX;
  do {
    // out += f(i)*(x - i)^(p - 1)
    deltaX = PowerXPlusAMod(-index, modulus - 1, polyModulus);
    deltaX = 1 - deltaX;
    out += f(index) * deltaX;

    index += 1;
  } while (!IsZero(index));

  return out;
  
}

NTL::ZZ_p equalTo(NTL::ZZ_p input) {
  return NTL::ZZ_p(input == NTL::ZZ_p(5));
}

Yashe::Ciphertext_t evaluatePoly(
    Yashe & SHE,
    NTL::ZZ_pX poly, 
    Yashe::Ciphertext_t input) {

  long sqrtDegree = sqrt(deg(poly));

  while ( (sqrtDegree + 1) * sqrtDegree < deg(poly) + 1 ) {
    sqrtDegree += 1;
  }

  std::vector<Yashe::Ciphertext_t> powers(sqrtDegree + 1);
  powers[0] = NTL::ZZ_pX(1);
  powers[1] = input;
  for (long i = 2; i <= sqrtDegree; i++) {
    if (i % 2 == 0) {
      powers[i] = SHE.mul(powers[i/2], powers[i/2]);
    } else {
      powers[i] = SHE.mul(powers[i - 1], input);
    }
  }


  std::vector<Yashe::Ciphertext_t> powersOfPowers(sqrtDegree + 1);
  powersOfPowers[0] = NTL::ZZ_p(1);
  powersOfPowers[1] = powers[sqrtDegree];
  for (long i = 2; i <= sqrtDegree; i++) {
    if (i % 2 == 0) {
      powersOfPowers[i] = SHE.mul(powersOfPowers[i/2], powersOfPowers[i/2]);
    } else {
      powersOfPowers[i] = SHE.mul(powersOfPowers[i - 1], powers[sqrtDegree]);
    }
  }

  Yashe::Ciphertext_t total = NTL::ZZ_pX(0);
  for (long chunk = 0; chunk <= sqrtDegree; chunk++) {
    Yashe::Ciphertext_t subtotal = NTL::ZZ_pX(0);
    for (long i = 0; 
         i < std::min(sqrtDegree, deg(poly) + 1 - chunk * sqrtDegree); 
         i++) {
      subtotal += SHE.mulcst(powers[i], rep(poly[i + chunk * sqrtDegree]) % deg(poly) + 1);
    }
    if (chunk > 0) {
      total += SHE.mul(subtotal, powersOfPowers[chunk]);
    } else {
      total = subtotal;
    }
  }

  return total;

}

NTL::ZZ_p evaluatePoly(NTL::ZZ_pX poly, NTL::ZZ_p input) {

  long sqrtDegree = sqrt(deg(poly));

  while ( (sqrtDegree + 1) * sqrtDegree < deg(poly) + 1 ) {
    sqrtDegree += 1;
  }

  std::vector<NTL::ZZ_p> powers(sqrtDegree + 1);
  powers[0] = NTL::ZZ_p(1);
  for (long i = 1; i <= sqrtDegree; i++) {
    if (i % 2 == 0) {
      powers[i] = powers[i/2] * powers[i/2];
    } else {
      powers[i] = powers[i - 1] * input;
    }
  }


  std::vector<NTL::ZZ_p> powersOfPowers(sqrtDegree + 1);
  powersOfPowers[0] = NTL::ZZ_p(1);
  for (long i = 1; i <= sqrtDegree; i++) {
    if (i % 2 == 0) {
      powersOfPowers[i] = powersOfPowers[i/2] * powersOfPowers[i/2];
    } else {
      powersOfPowers[i] = powersOfPowers[i - 1] * powers[sqrtDegree];
    }
  }

  NTL::ZZ_p total(0);
  for (long chunk = 0; chunk <= sqrtDegree; chunk++) {
    NTL::ZZ_p subtotal(0);
    for (long i = 0; 
         i < std::min(sqrtDegree, deg(poly) + 1 - chunk * sqrtDegree); 
         i++) {
      subtotal += poly[i + chunk * sqrtDegree] * powers[i];
    }
    total += subtotal * powersOfPowers[chunk];
  }

  return total;

}

int main() {

  // Cipher text modulus
  NTL::ZZ q = NTL::power2_ZZ(521) - 1;
  // Plain text modulus
  NTL::ZZ t = NTL::conv<NTL::ZZ>(12);
  // ??
  uint64_t d = 8;
  // Standard deviation
  uint64_t sigma = 8;
  // ??
  uint64_t base = 2; // 2^32
  
  Yashe SHE(q, t, d, sigma, base);

  SHE.keygen();

  std::cout << "here" << std::endl;

  //test_add(SHE);

  Yashe::Plainetxt_t m = NTL::ZZ_pX(5);
  std::cout << "message: " <<m << std::endl;

  Yashe::Ciphertext_t c = SHE.encrypt(m);
  std::cout << "ciphertext: "<< c << std::endl;

  Yashe::Plainetxt_t dec = SHE.decrypt(c);
  std::cout << "decryption: " << dec << std::endl;




  //test_add(SHE);
  //test_mul(SHE);
  //test_add_int64(SHE);

  //NTL::ZZ modulus = NTL::ZZ(7);

  //NTL::ZZ_p::init(modulus);

  //NTL::ZZ_pX poly = functionToCoefficients(equalTo);

  //std::cout << poly << std::endl;

  //for (long i = 0; i < 19; i++) {
    //std::cout << evaluatePoly(poly, NTL::ZZ_p(i)) << std::endl;
  //}

  return 0;
}
