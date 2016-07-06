#include <iostream>
#include <random>
#include <vector>
#include <functional>

#include "../source/functions.hpp"

bool testDivisionPolynomial() {
  
  srand(time(0));

  long t = 257;

  long denominator = rand() % t;
  if (denominator == 0) {
    denominator += 1;
  }

  std::function<long(long)> divideByDenom = Functions::divideByConstant(denominator);

  std::vector<long> poly = Functions::functionToPoly(divideByDenom, t);

  long numerator = rand() % t;

  long result = 0;
  for (long i = poly.size() - 1; i >= 0; i--) {
    result = (result * numerator + poly[i] ) % t;
  }

  return result == numerator/denominator;

}
  

int main() {
  if (not testDivisionPolynomial() ) {
    std::cout << "Testing division polynomial failed!" << std::endl;
  } else {
    std::cout << "Passed: division polynomial" << std::endl;
  }
  return 0;
}
