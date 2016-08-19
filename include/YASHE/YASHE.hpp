#ifndef YASHE_H
#define YASHE_H

#include <vector>
#include <random>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class YASHE_CT;

class YASHE {
  private:
    // Some stuff to write YASHE to file
    // using boost serialization.
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
        ar & pModulus;
        ar & bigPModulus;
        ar & cModulus;
        NTL::ZZ_p::init(cModulus);
        ar & bigModulus;
        ar & modulusRatio;
        ar & cycloModX;
        ar & maxDegree;
        ar & factors;
        ar & radix;
        ar & decompSize;
      }

    // The plain text modulus t
    long pModulus;
    // The plain text modulus t
    // wrapped as an NTL::ZZ int
    NTL::ZZ bigPModulus;
    // The cipher text modulus q
    NTL::ZZ cModulus;
    // The modulus used for multiplication, q*q/t
    NTL::ZZ bigModulus;
    // The ratio of the plain text
    // and cipher text moduli, q/t
    NTL::ZZ_p modulusRatio;

    // The dth cyclotomic polynomial
    // Phi_d(X)
    NTL::ZZX cycloModX;
    // The dth cyclotomic polynomial
    // modulo q
    NTL::ZZ_pXModulus cycloMod;
    // The dth cyclotomic polynomial
    // modulo q*q/t
    NTL::ZZ_pXModulus bigCycloMod;
    // The maximum degree of any polynomial
    // modulo the dth cyclotomic polynomial
    // = phi(d) - 1
    long maxDegree;

    // The factors of the dth cyclotomic polynomial
    std::vector<NTL::ZZ_pX> factors;

    // The radix w used to compute powers of Radix
    // or Radix decomp
    NTL::ZZ radix;
    // The size of decompositions, l_w,q = log_w(q) + 1
    long decompSize;


    // The standard deviation of the Gaussian
    // random variable, sigma.
    long stdDev;
    // A random number generator
    std::mt19937 randGen;

    // The public key
    NTL::ZZ_pX publicKey;
    // The evaluation key, with pro computation in order
    // to perform fast multiplication.
    std::vector<NTL::ZZ_pXMultiplier> evalKeyMult;

    std::vector<NTL::ZZ_pX> crtElements;
  public:
    /**
     * A constructor for the YASHE class
     * that requires all necessary parameters,
     * the plain text modulus, cipher text modulus,
     * cyclotomic index, standard deviation,
     * and decomposition base.
     *
     * These parameters will depend on the security 
     * parameter lambda, as well as the depth of computation
     * desired. q must be prime. 
     *
     * For a multiplicative depth of about 10, 8 bits of 
     * plain text space, a batch size of 10752, and a security
     * parameter of about 80, I recommend:
     *
     *    t = 257
     *    q = a prime with 1600 bits
     *    d = 66048 = 2^9*3*43 - optimal for batch size
     *    sigma = 8
     *    radix = 2^300
     */
    YASHE(
          long pModulus_,    // Plain text modulus t
          long log2CModulus_, // Cipher text modulus q
          long cyclotomicDegree_, // the index of the cyclotomic polynomial
          long stdDev_,      // Standard deviation of the gaussian distribution
          long log2Radix_        // base of decomposition
        );

    YASHE();

    static std::string getFileName(
          long pModulus_,    // Plain text modulus t
          long log2CModulus_, // Cipher text modulus q
          long cyclotomicDegree_, // the index of the cyclotomic polynomial
          long stdDev_,      // Standard deviation of the gaussian distribution
          long log2Radix_        // base of decomposition
        );

    /**
     * A constructor that reads the class
     * from a file. This means a class can be
     * reused without having to regenerate
     * factors or primes.
     * The random number generator is reseeded.
     */
    static YASHE readFromFile(std::string filename);

    /**
     * Writes a class to a file.
     */
    void writeToFile(std::string filename);

    /**
     * A variety of public functions
     * to fetch private variables
     */
    long getPModulus();
    NTL::ZZ getCModulus();
    NTL::ZZ getBigModulus();
    NTL::ZZ_p getModulusRatio();
    long getNumFactors();
    long getMaxDegree();
    NTL::ZZX getCycloModX();
    NTL::ZZ_pXModulus getCycloMod();
    NTL::ZZ_pXModulus getBigCycloMod();


    /**
     * returns a random polynomial with coefficients
     * in {-1,0,1}
     */
    NTL::ZZ_pX randomKeyPoly();

    /**
     * Returns a random polynomial with coefficients
     * as sampled from a Gaussian with standard
     * deviation sigma, centered around 0
     */
    NTL::ZZ_pX randomErrPoly();

    /**
     * Computes the dot product of two
     * vectors of polynomials
     */
    void dot(NTL::ZZ_pX& output,
             const std::vector<NTL::ZZ_pX>& a,
             const std::vector<NTL::ZZ_pX>& b);

    /**
     * Computes the dot of a with the evaluation key
     * This is bulk of the key switching step and
     * is optimized with pro computation on the eval
     * key.
     */
    void dotEval(NTL::ZZ_pX& output,
             const std::vector<NTL::ZZ_pX>& a);

    /**
     * Decomposes an integer into its digits base w:
     *
     *   input = output[0] + ouput[1]*w + output[2]*w^2 ...
     *
     */
    void radixDecomp(std::vector<NTL::ZZ> &output,
                     const NTL::ZZ& integer);
    /**
     * Decomposes a polynomial into a series of polynomial
     * base w:
     *
     *   input = output[0] + ouput[1]*w + output[2]*w^2 ...
     *
     */
    void radixDecomp(std::vector<NTL::ZZ_pX>& output,
                     const NTL::ZZ_pX& poly);

    /**
     * Returns a vector of the input multiplied by 
     * sucessive powers of w:
     *
     *   output = {input, input * w, input * w^2 ...}
     *
     * Note that the dot product of powersOfRadix(x) and
     * radixDecomp(x) is the original input x.
     */
    void powersOfRadix(std::vector<NTL::ZZ_pX> & output,
                       const NTL::ZZ_pX & poly);

    /**
     * Performs the rounded multiplication done in
     * the multiplication step:
     *
     *     ~c_mult = round(t/q*a*b) mod q
     * 
     * Note that the product within the round is not
     * taken to be modulo q.
     */
    void roundMultiply(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);
    /**
     * Performs round multiplication but when the
     * right input has already been preprocessed.
     * This is faster than regular round multiplication
     */
    void roundMultiply(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pXMultiplier& b);
    void roundMultiply(NTL::ZZ_pX& output, const NTL::ZZ_pXMultiplier& a, const NTL::ZZ_pXMultiplier& b);

    /**
     * Performs the round decryption step
     * that happens during decryption:
     *
     *    m = round(t/q*(a*b mod q)) mod t
     *
     */
    void roundDecrypt(NTL::ZZ& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);
    /**
     * Performs the round decryption step
     * on an entire vector, not just the first
     * integer.
     */
    void roundDecryptVec(NTL::ZZ_pX& output, const NTL::ZZ_pX& a, const NTL::ZZ_pX& b);

    /**
     * Performs a key switching operation on the input
     * and places the result in output.
     */
    void keySwitch(NTL::ZZ_pX& output, const NTL::ZZ_pX& input);

    /**
     * Generates public, secret and evaluation keys.
     * The secret key is returned to be stored by the user,
     * the public and evaluation keys are kept internally.
     */
    NTL::ZZ_pX keyGen();

    /**
     * Encrypts an integer message. The message
     * must be modulo t, the plain text modulus
     */
    YASHE_CT encrypt(long message);
    /**
     * Encrypts a vector of integer messages. The
     * messages must be modulo t, the plain text modulus.
     * All operations done to the ciphertext operate
     * in parallel on the vector of messages.
     * The input vector must be of size factors.length()
     */
    YASHE_CT encryptBatch(std::vector<long> messages);

    /**
     * Decrypts a single integer message
     */
    long decrypt(YASHE_CT ciphertext, NTL::ZZ_pX secretKey);
    /**
     * Decrypts and returns a vector of integer messages.
     * The resulting vector is of size factors.length()
     */
    std::vector<long> decryptBatch(YASHE_CT ciphertext, NTL::ZZ_pX secretKey);
};

#endif
