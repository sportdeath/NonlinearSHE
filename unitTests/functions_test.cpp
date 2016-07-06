#include <iostream>
#include <random>
#include <vector>
#include <functional>

#include "gtest/gtest.h"

#include "../source/functions.hpp"

TEST(PolynomialTest, DivisionPolynomial) {
  srand(time(0));

  long t = 257;

  // Do 10 random tests
  for (long i = 0; i < 10; i++) {
    long numerator = rand() % t;
    long denominator = rand() % t;
    if (denominator == 0) {
      denominator += 1;
    }

    // create division polynomial
    std::function<long(long)> divideByDenom = Functions::divideByConstant(denominator);
    std::vector<long> poly = Functions::functionToPoly(divideByDenom, t);

    // Evaluate the polynomial
    long result = 0;
    for (long i = poly.size() - 1; i >= 0; i--) {
      result = (result * numerator + poly[i] ) % t;
    }

    ASSERT_EQ(result, numerator/denominator);
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
