#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>
#include <functional>

class Functions {
  public:

    /**
     * Takes an input function on integers and turns it
     * into a polynomial that evaluates to that function
     * in the range 0... modulus - 1. 
     * Modulus must be prime.
     */
    static std::vector<long> functionToPoly(std::function<long(long)> f, long modulus);

    // division by constant
    static std::function<long(long)> divideByConstant(long denominator);

    // inequalities
    static std::function<long(long)> gt(long rightSide);
    static std::function<long(long)> lt(long rightSide);
    static std::function<long(long)> eq(long rightSide);
    static std::function<long(long)> leq(long rightSide);
    static std::function<long(long)> geq(long rightSide);
};


#endif
