#ifndef YASHE_H
#define YASHE_H

#include <vector>
#include <random>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

class YASHE_CT;

class YASHE {
  private:
    long pModulus;          // t
    NTL::ZZ bigPModulus;
    NTL::ZZ cModulus;       // q
    NTL::ZZ bigModulus;     // q*q/t
    NTL::ZZ_p modulusRatio;

    long cycloDegree;   // d
    long maxDegree;         // phi(d) - 1
    NTL::ZZX cycloModX;         // Phi_d(X)
    NTL::ZZ_pXModulus cycloMod; // Phi_d(X)
    NTL::ZZ_pXModulus bigCycloMod; // Phi_d(X)

    NTL::vec_ZZ_pX factors;

    NTL::ZZ radix;             // w
    long decompSize;        // l_w,q

    long stdDev;            // sigma
    std::mt19937 randGen;   // random number generator

    NTL::ZZ_pX publicKey;
    std::vector<NTL::ZZ_pXMultiplier> evalKeyMult;

  public:
    YASHE(
          long pModulus_,    // Plain text modulus t
          NTL::ZZ cModulus_, // Cipher text modulus t
          long degree_,      // degree of polynomial
          long stdDev_,      // Standard deviation of noise
          NTL::ZZ radix_        // base of decomp w
        );

    long getPModulus();
    NTL::ZZ getCModulus();
    NTL::ZZ getBigModulus();
    NTL::ZZ_p getModulusRatio();
    long getNumFactors();
    long getMaxDegree();
    NTL::ZZ_pXModulus getCycloMod();
    NTL::ZZ_pXModulus getBigCycloMod();

    NTL::ZZ_pX randomKeyPoly();

    NTL::ZZ_pX randomErrPoly();

    void dot(NTL::ZZ_pX& output,
             const std::vector<NTL::ZZ_pX>& a,
             const std::vector<NTL::ZZ_pX>& b);

    void dotEval(NTL::ZZ_pX& output,
             const std::vector<NTL::ZZ_pX>& a);

    void radixDecomp(std::vector<NTL::ZZ> &output,
                     const NTL::ZZ& integer);
    void radixDecomp(std::vector<NTL::ZZ_pX>& output,
                     const NTL::ZZ_pX& poly);

    void powersOfRadix(std::vector<NTL::ZZ_pX> & output,
                       const NTL::ZZ_pX & poly);

    void roundMultiply(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);
    void roundMultiply(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pXMultiplier& b);
    void roundDecrypt(NTL::ZZ& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);
    void roundDecryptVec(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);

    void keySwitch(NTL::ZZ_pX& output, const NTL::ZZ_pX& input);

    /**
     * Generates public, secret and evaluation keys.
     * The secret key is returned to be stored by the user.
     */
    NTL::ZZ_pX keyGen();

    YASHE_CT encrypt(std::vector<long> message);
    YASHE_CT encrypt(long message);

    YASHE_CT encryptBatch(std::vector<long> messages);

    std::vector<long> decryptVec(YASHE_CT ciphertext, NTL::ZZ_pX secretKey);
    long decrypt(YASHE_CT ciphertext, NTL::ZZ_pX secretKey);
    std::vector<long> decryptBatch(YASHE_CT ciphertext, NTL::ZZ_pX secretKey);
};

#endif
