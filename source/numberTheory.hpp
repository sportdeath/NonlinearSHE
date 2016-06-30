#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H

#include <NTL/ZZ_pX.h>

class NumberTheory {
  public:
    static long GCD(long a, long b);

    static long eulerToitient(long n);

    static NTL::ZZX cyclotomicPoly(long n);

    static void CRT(
        NTL::ZZ_pX& output,
        const std::vector<long>& inputs,
        const NTL::ZZ_pX& modulus,
        const NTL::vec_ZZ_pX& factors);
};

#endif
