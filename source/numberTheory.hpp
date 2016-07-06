#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H

#include <NTL/ZZ_pX.h>

/**
 * A class of basic number theoretic
 * functions which are used by the YASHE
 * scheme.
 */
class NumberTheory {
  public:
    /**
     * Returns the GCD of a and b
     */
    static long GCD(long a, long b);

    /**
     * Returns Euler's Totient function
     * evaluated at n. The value of Euler's
     * Totient function is equal to the number
     * of positive integers less than n
     * that are relatively prime to n.
     */
    static long eulerToitient(long n);

    /**
     * Returns the nth cyclotomic polynomial.
     * The nth cyclotomic polynomial is the
     * unique irriducible polynomial with integer
     * coefficients, which is a divisor of
     * x^n - 1 and not of x^k - 1 for all k < n.
     *
     * The degree of the nth cyclotomic polynomial
     * is equal to the Euler Totient function
     * evaluated at n.
     */
    static NTL::ZZX cyclotomicPoly(long n);

    /**
     * This computes a polynomial whose existence
     * is given by the Chinese Remainder theorem.
     * This polynomial (output) is the unique polynomial
     * that obeys the following condition
     *
     *        output = inputs[i] modulo factors[i]
     * 
     * for all i. This is very useful for batch encryption
     * of ciphertexts. The modulus is the product of all
     * of the factors, given in order to reduce running
     * time of computation.
     */
    static void CRT(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const std::vector<NTL::ZZ_pX> & factors);

    static void CRTwithElements(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const std::vector<NTL::ZZ_pX> & crtElements);
};

#endif
