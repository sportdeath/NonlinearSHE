#ifndef CIPHERTEXT_H
#define CIPHERTEXT_H

#include <vector>

#include <NTL/ZZ_pX.h>

#include "yashe.hpp"

class YASHE_CT {
  private:
    // A reference to the encryption scheme
    YASHE * y;

    // The polynomial representing the ciphertext
    NTL::ZZ_pX poly;

    // A precomputed value to make multiplication faster
    NTL::ZZ_pXMultiplier multiplier;

    // True if and only if the multiplier has been set
    bool isMultiplier = false;

  public:
    /**
     * Creates a blank ciphertext encrypting 0
     * Note that this ciphertext is not secure as
     * there is no added noise.
     */
    YASHE_CT();
    /**
     * Creates a ciphertext encrypting value "value"
     * under encryption scheme y_. Note that this
     * ciphertext is not secure as there is no
     * added noise.
     */
    YASHE_CT(long value, YASHE * y_);
    /**
     * Creates a ciphertext with ciphertext poly
     * under encryption scheme y_. Note that poly
     * is taken to be the ciphertext itself and
     * no further encryption occurs. To encrypt, use
     * the encrypt or encryptBatch functions of the
     * YASHE class.
     */
    YASHE_CT(NTL::ZZ_pX poly_, YASHE * y_);

    /**
     * Returns the polynomial representing the
     * ciphertext.
     */
    NTL::ZZ_pX getPoly();

    /**
     * If this particular ciphertext is going multiplied lots of times
     * (as is the case with the input to a polynomial...)
     * pre computation (Fourier transform) is done to make to
     * make multiplication faster.
     */
    void generateMultiplier();
    
    /**************************************
     * Procedural forms of basic arithmetic
     **************************************/
    static void add(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);
    static void sub(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);
    static void mul(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);

    static void mul(YASHE_CT& output, const YASHE_CT& a, const long& b);
    static void mul(YASHE_CT& output, const long& a, const YASHE_CT& b);

    static void add(YASHE_CT& output, const YASHE_CT& a, const long& b);
    static void add(YASHE_CT& output, const long& a, const YASHE_CT& b);

    /**
     * Evaluates the polynomial "poly" on input
     * "input" and returns the value to "output".
     * The polynomial is represented as a series of
     * coefficients. poly[i] is is the ith coefficient.
     *
     * output = poly[0] + input * poly[1] + input^2 * poly[2] ...
     */
    static void evalPoly(YASHE_CT& output,
                         YASHE_CT& input,
                         const std::vector<long>& poly
                         );

    /**
     * Computes homomorphic division
     */
    static void div(YASHE_CT& output, YASHE_CT& a, YASHE_CT& b);
};

#endif
