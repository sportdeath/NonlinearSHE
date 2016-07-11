#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>
#include <functional>

#include <NTL/ZZX.h>

class Functions {
  public:

    /**
     * Takes an input function on integers and turns it
     * into a polynomial that evaluates to that function
     * in the range 0... modulus - 1. The polynomial will
     * be of degree modulus - 1.
     * The modulus must be prime.
     */
    static NTL::ZZ_pX functionToPoly(std::function<long(long)> f, long modulus);

    /*********************************
     * These are some basic functions that
     * can be transformed into polynomials
     * for homomorphic operation
     *********************************/

    // Division by a constant
    static std::function<long(long)> divideByConstant(long denominator);

    // Inequalities.
    static std::function<long(long)> gt(long rightSide);
    static std::function<long(long)> lt(long rightSide);
    static std::function<long(long)> eq(long rightSide);
    static std::function<long(long)> leq(long rightSide);
    static std::function<long(long)> geq(long rightSide);
};


#endif
