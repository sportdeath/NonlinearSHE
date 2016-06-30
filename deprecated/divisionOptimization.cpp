void YASHE_CT::div(YASHE_CT& output, const YASHE_CT& a, const YASHE_CT& b) {
  // computes the log base 2 of the input.
  // Only need log_2log_2 p integer bits
  // log_2 p - log_2log_2 p - 1 floating bits
  // 1 bit which will become high if the input is less than the output

  long t = a.y -> getPModulus();

  double base = 1.9;
  long minimumCount = 256*257;
  long miniMaxError = 256;
  double minimumBase = 1.9;
  double minimumBaseE = 1.9;
  for (long i = 0; i < 1000; i++) {
    base += 0.001;
  std::function<long(long)> divLog2 = [t, base](long input) {

    if (input == 0) {
      return long(0);
    }

    double logInput = log(input)/log(base);

    double shiftAmmount = log(t)/log(base) - log(log(t)/log(base))/log(base) - 1;

    long output = round(logInput * pow(base, shiftAmmount));

    return output;
  };


  std::function<long(long)> divExp2 = [t, base](long input) {
    double maxInt = log(t)/log(base);
    double shiftAmmount = log(t)/log(base) - log(log(t)/log(base))/log(base) - 1;
    double actualValue = input / pow(base, shiftAmmount);

    if (actualValue > maxInt) {
      return long(0);
    } else {
      return long(pow(base, actualValue));
    }
  };
    long maxError = 0;
    long maxN, maxD;
    long count = 0;
    for (long numerator = 0; numerator < t; numerator ++) {
      for (long denominator = 1; denominator < t; denominator ++) {
        long exponent = (divLog2(numerator) - divLog2(denominator)) % t;
        while (exponent < 0) exponent += t;
        long result = divExp2(exponent);

        if (result != numerator/denominator) {
          count += 1;
          long error = std::abs(numerator/denominator - result);
          if (error > maxError) {
            maxError = error;
            maxN = numerator;
            maxD = denominator;

          }
        }
      }
    }
    std::cout << count/(256. * 257.)*100. << "% incorrect" << std::endl;
    std::cout << "max error: " << maxError << " with " << maxN <<"/" << maxD << std::endl;
    std::cout << "base: " << base << std::endl;
    if (count < minimumCount) {
      minimumCount = count;
      minimumBase = base;
    }
    if (maxError < miniMaxError) {
      miniMaxError = maxError;
      minimumBaseE = base;
    }
  }
  std::cout << "minimum: " << minimumCount/(256.*257.)*100. << "% incorrect" << std::endl;
  std::cout << "with base: " << minimumBase << std::endl;
  std::cout << "minimum: " << miniMaxError << std::endl;
  std::cout << "with base: " << minimumBaseE << std::endl;

}
