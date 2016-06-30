#ifndef CIPHERTEXT_H
#define CIPHERTEXT_H

#include <vector>

#include <NTL/ZZ_pX.h>

#include "yashe.hpp"

class YASHE_CT {
  private:
    YASHE * y;
    NTL::ZZ_pX poly;
    NTL::ZZ_pXMultiplier multiplier;
    bool isMultiplier = false;

  public:
    YASHE_CT();
    YASHE_CT(long value, YASHE * y_);
    YASHE_CT(NTL::ZZ_pX poly_, YASHE * y_);

    NTL::ZZ_pX getPoly();

    /**
     * If this particular is going multiplied lots of times
     * (as is the case with the input to a polynomial...)
     * pre computation (Fourier transform) is done to make to
     * make multiplication faster
     */
    void generateMultiplier();
    
    // Procedural forms:
    static void add(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);
    static void sub(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);
    static void mul(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b);

    static void mul(YASHE_CT& output, const YASHE_CT& a, const long& b);
    static void mul(YASHE_CT& output, const long& a, const YASHE_CT& b);

    static void add(YASHE_CT& output, const YASHE_CT& a, const long& b);
    static void add(YASHE_CT& output, const long& a, const YASHE_CT& b);

    static void evalPoly(YASHE_CT& output,
                         YASHE_CT& input,
                         const std::vector<long>& poly
                         );

    static void div(YASHE_CT& output, YASHE_CT& a, YASHE_CT& b);
};

#endif
