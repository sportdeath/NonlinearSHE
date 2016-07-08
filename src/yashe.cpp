#include <vector>
#include <random>
#include <algorithm>
#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>

#include<boost/archive/text_oarchive.hpp>
#include<boost/archive/text_iarchive.hpp>
#include<boost/serialization/utility.hpp>
#include<boost/serialization/vector.hpp>

#include "YASHE/serialization/ZZ.hpp"
#include "YASHE/serialization/ZZX.hpp"
#include "YASHE/serialization/ZZ_p.hpp"
#include "YASHE/serialization/ZZ_pX.hpp"
#include "YASHE/serialization/ZZ_pXModulus.hpp"
#include "YASHE/serialization/FFTRep.hpp"

#include "YASHE/YASHE.hpp"
#include "YASHE/cipherText.hpp"
#include "YASHE/numberTheory.hpp"

YASHE::YASHE() {
  randGen = std::mt19937(time(0));
}


YASHE::YASHE(long pModulus_,
             NTL::ZZ cModulus_,
             long cyclotomicDegree,
             long stdDev_,
             NTL::ZZ radix_) {
  pModulus = pModulus_;
  cModulus = cModulus_;
  bigPModulus = NTL::ZZ(pModulus);

  NTL::ZZ_p::init(cModulus);
  modulusRatio = NTL::conv<NTL::ZZ_p>(cModulus/pModulus);
  bigModulus = (cModulus * cModulus)/pModulus;

  maxDegree = NumberTheory::eulerToitient(cyclotomicDegree) - 1;
  cycloModX = NumberTheory::cyclotomicPoly(cyclotomicDegree);
  cycloMod = NTL::ZZ_pXModulus(NTL::conv<NTL::ZZ_pX>(cycloModX));

  {
    NTL::ZZ_pPush push(bigPModulus); // switch to plain text modulus
    // Factor the cyclotomic polynomial modulo t
    // for batch encryption
    NTL::ZZ_pXModulus pModulusX;
    NTL::build(pModulusX, NTL::conv<NTL::ZZ_pX>(cycloModX));
    NTL::vec_ZZ_pX factors_;
    NTL::SFCanZass(factors_, pModulusX , 1);
    factors.resize(factors_.length());
    for (long i = 0; i < factors_.length(); i++) {
      factors[i] = factors_[i];
    }

    //NTL::ZZ_pX fInv, fInvInv;
    //crtElements.resize(factors.size());

    //for (long i = 0; i < factors.size(); i++) {
      //div(fInv, pModulusX, factors[i]);
      //rem(fInvInv, fInv, factors[i]);
      //InvMod(fInvInv, fInvInv, factors[i]);

      //crtElements[i] = MulMod(fInv, fInvInv, pModulusX);
    //}
  }

  {
    NTL::ZZ_pPush push(bigModulus); // switch to multiplication modulus
    // make another modulus for fast multiplication
    bigCycloMod = NTL::ZZ_pXModulus(NTL::conv<NTL::ZZ_pX>(cycloModX));
  }

  stdDev = stdDev_;
  radix = radix_;
  decompSize = log(cModulus)/log(radix) + 1; // log_w(q) + 1

  randGen = std::mt19937(time(0));

}


void YASHE::writeToFile(std::string filename) {
  std::ofstream ofs(filename);
  boost::archive::text_oarchive oa(ofs);
  oa << *this;
}

YASHE YASHE::readFromFile(std::string filename) {
  YASHE output;
  std::ifstream ifs(filename);
  boost::archive::text_iarchive ia(ifs);
  ia >> output;
  NTL::ZZ_p::init(output.cModulus);
  output.cycloMod = NTL::ZZ_pXModulus(NTL::conv<NTL::ZZ_pX>(output.cycloModX));
  {
    NTL::ZZ_pPush push(output.bigModulus); // switch to multiplication modulus
    // make another modulus for fast multiplication
    output.bigCycloMod = NTL::ZZ_pXModulus(NTL::conv<NTL::ZZ_pX>(output.cycloModX));
  }
  return output;
}


long YASHE::getPModulus() {
  return pModulus;
}


long YASHE::getNumFactors() {
  return factors.size();
}


NTL::ZZ YASHE::getCModulus() {
  return cModulus;
}

NTL::ZZ YASHE::getBigModulus() {
  return bigModulus;
}

NTL::ZZ_p YASHE::getModulusRatio() {
  return modulusRatio;
}

long YASHE::getMaxDegree() {
  return maxDegree;
}

NTL::ZZ_pXModulus YASHE::getCycloMod() {
  return cycloMod;
}

NTL::ZZ_pXModulus YASHE::getBigCycloMod() {
  return bigCycloMod;
}

NTL::ZZ_pX YASHE::randomKeyPoly() {
  NTL::ZZ_pX output;
  output.SetLength(maxDegree + 1);

  for (long i = 0; i <= maxDegree; i++) {
    long choice = (randGen() % 3);
    choice -= 1; // choice is in {-1, 0, 1}
    output[i] = NTL::ZZ_p(choice);
  }

  return output;
}


NTL::ZZ_pX YASHE::randomErrPoly() {

  std::normal_distribution<double> dist(0,stdDev);

  NTL::ZZ_pX output;
  output.SetLength(maxDegree + 1);

  for (long i = 0; i <= maxDegree; i++) {
    output[i] = NTL::ZZ_p(dist(randGen));
  }

  return output;
}


void YASHE::radixDecomp(std::vector<NTL::ZZ> & output,
                        const NTL::ZZ & input) {

  output.resize(decompSize);
  NTL::ZZ remainder;
  NTL::ZZ quotient, numerator;
  numerator = input;

  for (long i = 0; i < decompSize; i++) {
    DivRem(quotient, remainder, numerator, radix);
    output[i] = remainder;
    numerator = quotient;
  }
}


void YASHE::radixDecomp(std::vector<NTL::ZZ_pX> & output,
                        const NTL::ZZ_pX & poly) {

  output.resize(decompSize);

  for (long i = 0; i < decompSize; i++) {
    output[i].SetLength(NTL::deg(poly) + 1);
  }

  std::vector<NTL::ZZ> coefDecomp;
  for (long i = 0; i <= NTL::deg(poly); i++) {
    radixDecomp(coefDecomp, rep(poly[i]));
    for (long j = 0; j < decompSize; j++) {
      output[j][i] = NTL::conv<NTL::ZZ_p>(coefDecomp[j]);
    }
  }
}


void YASHE::powersOfRadix(std::vector<NTL::ZZ_pX> & output,
                          const NTL::ZZ_pX & poly) {

  output.resize(decompSize);

  output[0] = poly;

  for (long i = 1; i < decompSize; i++) {
    output[i] = output[i - 1] * NTL::conv<NTL::ZZ_pX>(radix);
  }
}


void YASHE::dot(NTL::ZZ_pX& output, 
                const std::vector<NTL::ZZ_pX>& a, 
                const std::vector<NTL::ZZ_pX>& b) {
  output = NTL::ZZ_pX::zero();
  NTL::ZZ_pX product;
  for (long i = 0; i < decompSize; i++) {
    NTL::MulMod(product, a[i], b[i], cycloMod);
    output += product;
  }
}


/**
 * Precomputation on the evaluation key
 * allows for faster key switching
 */
void YASHE::dotEval(NTL::ZZ_pX& output, 
                    const std::vector<NTL::ZZ_pX>& a) {
  long n = cycloMod.n;
  long da = maxDegree;

  NTL::FFTRep fftOutput(NTL::INIT_SIZE, cycloMod.k);

  ToFFTRep(fftOutput, NTL::ZZ_pX(0), cycloMod.k);

  NTL::ZZ_pX P1(NTL::INIT_SIZE, n), P2(NTL::INIT_SIZE, n);
  NTL::FFTRep R1(NTL::INIT_SIZE, cycloMod.l), R2(NTL::INIT_SIZE, cycloMod.l);

  for (long i = 0; i < decompSize; i++) {
    ToFFTRep(R1, a[i], cycloMod.l);
    mul(R2, R1, evalKeyMult[i].B1);
    FromFFTRep(P1, R2, n-1, 2*n-3);
    reduce(R1, R1, cycloMod.k);
    mul(R1, R1, evalKeyMult[i].B2);
    ToFFTRep(R2, P1, cycloMod.k);
    mul(R2, R2, cycloMod.FRep);
    sub(R1, R1, R2);

    add(fftOutput, R1, fftOutput);
  }

  FromFFTRep(output, fftOutput, 0, n-1);
}


void YASHE::roundMultiply(NTL::ZZ_pX& output,
                          const NTL::ZZ_pX& a,
                          const NTL::ZZ_pX& b) {

  // maximum is q^2 * (maxDegree + 1)
  NTL::ZZ_pPush push((cModulus * cModulus)/pModulus);


  NTL::ZZ_pX product;
  NTL::MulMod(product, a, b, bigCycloMod);

  output.SetLength(maxDegree + 1);
  NTL::ZZ quotient, remainder;

  for (long i = 0; i <= maxDegree; i++) {
    DivRem(quotient, remainder, pModulus * rep(product[i]), cModulus);

    // Rounding using remainder
    if (remainder * 2 > cModulus) {
      quotient += 1;
    }

    output[i] = NTL::conv<NTL::ZZ_p>(quotient);
  }
}


void YASHE::roundMultiply(NTL::ZZ_pX& output,
                          const NTL::ZZ_pX& a,
                          const NTL::ZZ_pXMultiplier& b) {
  // maximum is q^2 * (maxDegree + 1)
  NTL::ZZ_pPush push((cModulus * cModulus)/pModulus);

  NTL::ZZ_pX product;
  NTL::MulMod(product, a, b, bigCycloMod);

  output.SetLength(maxDegree + 1);
  NTL::ZZ quotient, remainder;

  for (long i = 0; i <= maxDegree; i++) {
    DivRem(quotient, remainder, pModulus * rep(product[i]), cModulus);

    // Rounding using remainder
    if (remainder * 2 > cModulus) {
      quotient += 1;
    }

    output[i] = NTL::conv<NTL::ZZ_p>(quotient);
  }
}

void YASHE::roundDecryptVec(NTL::ZZ_pX& output,
                            const NTL::ZZ_pX& a,
                            const NTL::ZZ_pX& b) {
  NTL::ZZ_pX product = MulMod(a, b, cycloMod);

  output.SetLength(maxDegree + 1);
  NTL::ZZ quotient, remainder;

  for (long i = 0; i <= maxDegree; i++) {

    DivRem(quotient, remainder, pModulus * rep(product[i]), cModulus);

    // Rounding using remainder
    if (remainder * 2 > cModulus) {
      quotient += 1;
    }

    output[i] = NTL::conv<NTL::ZZ_p>(quotient);
  }
}


void YASHE::roundDecrypt(NTL::ZZ& output,
                         const NTL::ZZ_pX& a,
                         const NTL::ZZ_pX& b) {
  NTL::ZZ_pX product = MulMod(a, b, cycloMod);

  NTL::ZZ quotient, remainder;

  DivRem(output, remainder, pModulus * rep(product[0]), cModulus);

    // Rounding using remainder
  if (remainder * 2 > cModulus) {
    output += 1;
  }
}


NTL::ZZ_pX YASHE::keyGen() {
  /**
   * The secret key is computed as
   *
   *      f' <- X_key
   *      f = (t*f' + 1) mod q
   *      secretKey = f
   * 
   * Secret keys are randomly generated until
   * an invertible key f^-1 is found
   */
  NTL::ZZ_pX secretKey, secretKeyInv;
  do {
    secretKey = pModulus * randomKeyPoly() + 1;
    secretKeyInv = PowerMod(secretKey, -1, cycloMod);
  } while (MulMod(secretKey, secretKeyInv, cycloMod) != NTL::ZZ_pX(1));

  /**
   * The public key is computed by
   *
   *      g <- X_key
   *      h = (t*g*f^-1) mod q
   *      publicKey = h
   */
  publicKey = MulMod(randomKeyPoly(), secretKeyInv, cycloMod);
  publicKey *= pModulus;

  /** 
   * The evaluation key is computed by
   *
   *        e, s <- X_err
   *        gamma = (powersOfRadix(f) + e + h*s) mod q
   *        evaluationKey = gamma
   */
  std::vector<NTL::ZZ_pX> evalKey;
  powersOfRadix(evalKey, secretKey);
  evalKeyMult.resize(decompSize);
  for (long i = 0; i < decompSize; i++) {
    evalKey[i] += randomErrPoly();
    evalKey[i] += MulMod(publicKey, randomErrPoly(), cycloMod);
    NTL::build(evalKeyMult[i], evalKey[i], cycloMod);
  }

  return secretKey;
}


/**
 * A message is a polynomial of degree <= d modulo t.
 * The encryption is calculated as follows:
 * 
 *   c = (floor(q/t) * (m mod t) + e + publicKey*s) mod q
 *
 * With s, e <- X_err
 */
YASHE_CT YASHE::encrypt(long message) {
  
  NTL::ZZ_pX output;
  output.SetLength(maxDegree + 1);

  output[0] = NTL::ZZ_p(message % pModulus);

  output *= modulusRatio;
  output += randomErrPoly();
  output += MulMod(randomErrPoly(), publicKey, cycloMod);

  return YASHE_CT(output, this);
}
YASHE_CT YASHE::encryptBatch(std::vector<long> messages) {
  
  NTL::ZZ_pX smallOutput;


  // Use the chinese remainder theorem to
  // encode the vector of messages into a
  // single polynomial.
  {
    NTL::ZZ_pPush push(bigPModulus);
    NumberTheory::CRT(smallOutput, messages, NTL::conv<NTL::ZZ_pX>(cycloModX), factors);
  }

  NTL::ZZ_pX output = NTL::conv<NTL::ZZ_pX>(NTL::conv<NTL::ZZX>(smallOutput));

  output *= modulusRatio;
  output += randomErrPoly();
  output += MulMod(randomErrPoly(), publicKey, cycloMod);

  return YASHE_CT(output, this);
}

/**
 * A cipher text is decrypted with the following
 * formula:
 *      
 *      m = round(t/q * ((f*c) mod q) ) mod t
 */
long YASHE::decrypt(YASHE_CT ciphertext, NTL::ZZ_pX secretKey) {
  
  NTL::ZZ decryption;
  roundDecrypt(decryption, secretKey, ciphertext.getPoly());

  return rem(decryption, pModulus);
}
std::vector<long> YASHE::decryptBatch(YASHE_CT ciphertext, NTL::ZZ_pX secretKey) {

  std::vector<long> output(factors.size());
  
  NTL::ZZ_pX decryption;
  roundDecryptVec(decryption, secretKey, ciphertext.getPoly());

  {
    NTL::ZZ_pPush push(bigPModulus);

    NTL::ZZ_pX smallDecryption;

    smallDecryption = NTL::conv<NTL::ZZ_pX>(decryption);

    // Undo the encoding of the chinese remainder theorem
    // to get each result out
    NTL::ZZ_pX remainder;
    for (long i = 0; i < factors.size(); i++) {
      rem(remainder, smallDecryption, factors[i]);
      output[i] = rem(rep(remainder[0]), pModulus);
    }
  }

  return output;
}


/**
 * Key switching is computed by
 *
 *     radixDecomp (dot) evalKey
 */
void YASHE::keySwitch(NTL::ZZ_pX& output, const NTL::ZZ_pX& input) {
  std::vector<NTL::ZZ_pX> decomp;
  radixDecomp(decomp, input);
  dotEval(output, decomp);
}
