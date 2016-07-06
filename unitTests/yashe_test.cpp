#include <vector>
#include <ctime>
#include <cmath>

#include <NTL/ZZ.h>

#include "../source/yashe.hpp"
#include "../source/cipherText.hpp"
#include "../source/functions.hpp"

#include "gtest/gtest.h"

class YASHETest : public ::testing::Test {
  protected:

    static long d, t;
    static NTL::ZZ q, w;
    static YASHE SHE;

};

long YASHETest::t = 257;
NTL::ZZ YASHETest::q = NTL::GenPrime_ZZ(400);
long YASHETest::d = 22016; // 2^9*43 - 5376 irreducible factors
NTL::ZZ YASHETest::w = NTL::power2_ZZ(70);
YASHE YASHETest::SHE = YASHE(t,q,d,8,w);

TEST_F(YASHETest, RandomKeyPolyBounds) {
  NTL::ZZ_pX randPoly = SHE.randomKeyPoly();
  for (long i = 0; i <= deg(randPoly) ) {
    ASSERT_LE(randPoly[i],1);
    ASSERT_GE(randPoly[i],-1);
  }
}

TEST_F(YASHETest, RandomKeyPolyDifferent) {
  NTL::ZZ_pX randPoly1 = SHE.randomKeyPoly();
  NTL::ZZ_pX randPoly2 = SHE.randomKeyPoly();

  ASSERT_NE(randPoly1, randPoly2);
}

TEST_F(YASHETest, RandomErrPolyBounds) {
  NTL::ZZ_pX randPoly = SHE.randomKeyPoly();

  for (long i = 0; i <= deg(randPoly) ) {
    ASSERT_LE(randPoly[i],1);
    ASSERT_GE(randPoly[i],-1);
  }
}

bool testRadixDecompInt() {
  YASHE SHE(2,NTL::ZZ(257),5,8,NTL::ZZ(3));

  std::vector<NTL::ZZ> decomp;
  SHE.radixDecomp(decomp, NTL::ZZ(200));

  
  std::vector<long> result {2, 0, 1, 1, 2, 0};

  bool isSame = true;

  for (long i = 0; i < decomp.size(); i++) {
    isSame = isSame & (decomp[i] == result[i]);
  }
  return isSame;
}

void testRadixDecompPoly() {
  YASHE SHE(2,NTL::ZZ(257),5,8,NTL::ZZ(3));

  NTL::ZZ_pX testPoly;
  testPoly.SetLength(5);
  std::vector<long> coefs {44, 126, 40, 131, 166};

  for (long i = 0; i < 5; i++) {
    testPoly[i] = coefs[i];
  }

  //                       44 126  40 131 166
  std::vector<long> out1 {  2,  0,  1,  2,  1};
  std::vector<long> out2 {  2,  0,  1,  1,  1};
  std::vector<long> out3 {  1,  2,  1,  2,  0};
  std::vector<long> out4 {  1,  1,  1,  1,  0};
  std::vector<long> out5 {  0,  1,  0,  1,  2};
  std::vector<long> out6 {  0,  0,  0,  0,  0};

  std::vector< std::vector<long> > outs = {out1, out2, out3, out4, out5, out6};

  std::vector<NTL::ZZ_pX> decomp;
  SHE.radixDecomp(decomp, testPoly);

  bool isSame = true;

  for (long i = 0; i < 6; i++) {
    for (long j = 0; j < 5; j++) {
      isSame = isSame & (decomp[i][j] == outs[i][j]);
    }
  }
      

  return isSame;
}

bool testPowersOfRadix() {
  YASHE SHE(2,NTL::ZZ(257),5,8,NTL::ZZ(3));

  NTL::ZZ_pX testPoly;
  testPoly.SetLength(5);
  std::vector<long> coefs {44, 126, 40, 131, 166};

  for (long i = 0; i < 5; i++) {
    testPoly[i] = coefs[i];
  }

  //                       44 126  40 131 166
  std::vector<long> out1 { 44,126, 40,131,166};
  std::vector<long> out2 {132,121,120,136,241}; // * 3
  std::vector<long> out3 {139,106,103,151,209}; // * 3^2
  std::vector<long> out4 {160, 61, 52,196,113}; // * 3^3
  std::vector<long> out5 {223,183,156, 74, 82}; // * 3^4
  std::vector<long> out6 {155, 35,211,222,246}; // * 3^5

  std::vector< std::vector<long> > outs = {out1, out2, out3, out4, out5, out6};

  std::vector<NTL::ZZ_pX> powers;
  SHE.powersOfRadix(powers, testPoly);

  bool isSame = true;

  for (long i = 0; i < 6; i++) {
    for (long j = 0; j < 5; j++) {
      isSame = isSame & (powers[i][j] == outs[i][j]);
    }
  }
  return isSame;
}

bool testDecompPowers() {

  YASHE SHE(2,NTL::ZZ(257),5,8,NTL::ZZ(3));

  NTL::ZZ_pX testPoly1 = SHE.randomKeyPoly();
  NTL::ZZ_pX testPoly2 = SHE.randomKeyPoly();

  NTL::ZZ_pX mulPoly = NTL::MulMod(testPoly1, testPoly2, SHE.getCycloMod());

  std::vector<NTL::ZZ_pX> powers, decomp;
  SHE.powersOfRadix(powers, testPoly1);
  SHE.radixDecomp(decomp, testPoly2);

  NTL::ZZ_pX dotPoly;
  SHE.dot(dotPoly, powers, decomp);

  return mulPoly==dotPoly;
}


bool testKeyGen() {
  YASHE SHE(2,NTL::ZZ(15688861),5,8,NTL::ZZ(3));

  SHE.keyGen();

  return true;
}

//bool testEncryptDecrypt() {
  //long t = 257;
  //NTL::ZZ q = NTL::GenPrime_ZZ(392);
  //long d = 16384;
  //NTL::ZZ w = NTL::power2_ZZ(32);
  //YASHE SHE(t,q,d,8,w);

  //NTL::ZZ_pX secretKey = SHE.keyGen();

  //std::vector<long> message(SHE.getMaxDegree() + 1);

  //srand(time(0));

  //for (long i = 0; i <= SHE.getMaxDegree(); i++) {
    //message[i] = rand() % t;
  //}

  //YASHE_CT ciphertext = SHE.encrypt(message);

  //std::vector<long> decryption = SHE.decryptVec(ciphertext, secretKey);

  //bool isSame = true;
  //for (long i = 0; i <=SHE.getMaxDegree(); i++) {
    //isSame &= (message[i] == decryption[i]);
  //}

  //return isSame;
//}


bool testEncryptDecryptBatch() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(392);
  long d = 2048;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  std::vector<long> message(SHE.getNumFactors());

  srand(time(0));

  for (long i = 0; i <= SHE.getNumFactors(); i++) {
    message[i] = rand() % t;
  }

  YASHE_CT ciphertext = SHE.encryptBatch(message);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

  bool isSame = true;
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    isSame &= (message[i] == decryption[i]);
  }

  return isSame;
}

bool testAddCiphertexts() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(392);
  long d = 16384;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  //std::vector<long> message1(SHE.getMaxDegree() + 1);
  //std::vector<long> message2(SHE.getMaxDegree() + 1);

  srand(time(0));

  long message1 = rand() % t;
  long message2 = rand() % t;

  //for (long i = 0; i <= SHE.getMaxDegree(); i++) {
    //message1[i] = rand() % t;
    //message2[i] = rand() % t;
  //}

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  clock_t begin = clock();

  YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);

  clock_t end = clock();

  std::cout << "addition took " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  //std::vector<long> decryption = SHE.decryptVec(ciphertext1, secretKey);
  long decryption = SHE.decrypt(ciphertext1, secretKey);

  //bool isSame = true;
  ////for (long i = 0; i <=SHE.getMaxDegree(); i++) {
    //isSame &= ((message1[i] + message2[i]) % t == decryption[i]);
  //}
  return (message1 + message2) == decryption;
}

bool testMulCiphertexts() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 32768;
  //long d = 256;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  long message1, message2;
  message1 = rand() % t;
  message2 = rand() % t;

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  clock_t begin = clock();

  YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);

  clock_t end = clock();

  std::cout << "multiplication took " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  long decryption = SHE.decrypt(ciphertext1, secretKey);


  return (message1 * message2) % t == decryption;
}


bool testAddBatchCiphertexts() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 2048;
  //long d = 256;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = rand() % t;
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  clock_t begin = clock();

  YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);

  clock_t end = clock();

  std::cout << "addition took " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  bool isSame = true;

  for (long i = 0; i < SHE.getNumFactors(); i++) {
    isSame &=  (message1[i] + message2[i]) % t == decryption[i];
  }

  return isSame;
}


bool testMulBatchCiphertexts() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 2048;
  //long d = 256;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = rand() % t;
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  clock_t begin = clock();

  YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);

  clock_t end = clock();

  std::cout << "multiplication took " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  bool isSame = true;

  for (long i = 0; i < SHE.getNumFactors(); i++) {
    isSame &=  (message1[i] * message2[i]) % t == decryption[i];
  }

  return isSame;
}


bool testAdditionByConstant() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 32768;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  long message = rand() % t;
  long constant = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);

  clock_t begin = clock();

  YASHE_CT::add(ciphertext, ciphertext, constant);

  clock_t end = clock();

  std::cout << "addition by constant took: " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  long decryption = SHE.decrypt(ciphertext, secretKey);

  return (message + constant) % t == decryption;
}


bool testMultiplicationByConstant() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 32768;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  long message = rand() % t;
  long constant = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);

  clock_t begin = clock();

  YASHE_CT::mul(ciphertext, ciphertext, constant);

  clock_t end = clock();

  std::cout << "multiplication by constant took: " << double(end - begin)/CLOCKS_PER_SEC << "seconds" << std::endl;

  long decryption = SHE.decrypt(ciphertext, secretKey);

  return (message * constant) % t == decryption;
}


bool testMultiplicativeDepth() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(800);
  long d = 32768;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  long depth = 0;

  long message1, message2, decryption;
  message1 = rand() % t;
  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  do {
    message2 = rand() % t;

    YASHE_CT ciphertext2 = SHE.encrypt(message2);

    YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);

    decryption = SHE.decrypt(ciphertext1, secretKey);

    message1 = (message1 * message2) % t;

    depth += 1;

    std::cout << message1 << " " << decryption << std::endl;
  } while (message1 == decryption);

  std::cout << "reached maximum multiplication depth of " << depth << "!" << std::endl;

  return true;
}

bool testAdditiveDepth() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(32);
  long d = 256;
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  long depth = 0;

  long message1, message2, decryption;
  message1 = rand() % t;
  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  do {
    message2 = rand() % t;

    YASHE_CT ciphertext2 = SHE.encrypt(message2);

    YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);

    decryption = SHE.decrypt(ciphertext1, secretKey);

    message1 = (message1 + message2) % t;

    depth += 1;

    std::cout << message1 << " " << decryption << std::endl;
  } while (message1 == decryption);

  std::cout << "reached maximum additive depth of " << depth << "!" << std::endl;

  return true;
}

bool testEvalPoly() {
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(1024);
  long d = 1024; 
  NTL::ZZ w = NTL::power2_ZZ(32);
  YASHE SHE(t,q,d,8,w);

  NTL::ZZ_pX secretKey = SHE.keyGen();

  srand(time(0));

  //std::cout << "poly = [";
  long degree = 256;
  std::vector<long> poly(degree);
  for (long i = 0; i < degree; i++) {
    poly[i] = rand() % t;
    //std::cout << poly[i] << " ";
  }
  //std::cout << "]" << std::endl;

  long message = rand() % t;

  //std::cout << "input: " << message << std::endl;

  YASHE_CT ciphertext = SHE.encrypt(message);

  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  long result = 0;
  for (long i = degree - 1; i >= 0; i--) {
    //if (i % 4 == 0) {
      //std::cout << "subtotal up to degree " << i - 1 << " = " << result << std::endl;
    //}
    result = (result * message + poly[i]) % t;
  }
  //std::cout << "desired result: " << result << std::endl;

  return result == decryption;

}

bool testDivisionByConstant() {


  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(1424);
  long d = 65536;
  NTL::ZZ w = NTL::power2_ZZ(300);
  YASHE SHE(t,q,d,8,w);

  clock_t start, end;

  start = clock();

  NTL::ZZ_pX secretKey = SHE.keyGen();

  end = clock();

  std::cout << "Keygen took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  srand(time(0));

  long denominator = rand() % t;
  if (denominator == 0) {
    denominator += 1;
  }

  std::vector<long> poly = Functions::functionToPoly(Functions::divideByConstant(denominator), t);

  long numerator = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(numerator);

  start = clock();

  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);
  end = clock();

  std::cout << "division took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  long decryption = SHE.decrypt(ciphertext, secretKey);

  return decryption == numerator/denominator;
}


bool testBatchDivisionByConstant() {
  long t = 257;
  //NTL::ZZ q = NTL::GenPrime_ZZ(400);
  ////long d = 688;
  //long d = 22016; // 2^9*43 - 5376 irreducible factors
  ////long d = 66048; // 2^9*3*43 - 10752 irreducible factors
  ////long d = 1376*2; // 2^5 * 43
  //NTL::ZZ w = NTL::power2_ZZ(70);
  //YASHE SHE(t,q,d,8,w);
  //SHE.writeToFile("8BitFHE");
  YASHE SHE = YASHE::readFromFile("/Users/tfh/Dropbox (MIT)/Work/Orange/YASHE/build/8BitFHE");


  std::cout << "batch size: " << SHE.getNumFactors() << std::endl;

  clock_t start, end;

  start = clock();

  NTL::ZZ_pX secretKey = SHE.keyGen();

  end = clock();

  std::cout << "Keygen took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  srand(time(0));

  long denominator = rand() % t;
  if (denominator == 0) {
    denominator += 1;
  }

  denominator = 2;

  std::vector<long> poly = Functions::functionToPoly(Functions::divideByConstant(denominator), t);

  std::vector<long> numerators(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    numerators[i] = rand() % t;
  }
  start = clock();

  YASHE_CT ciphertext = SHE.encryptBatch(numerators);
  end = clock();

  std::cout << "encryption took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  start = clock();

  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);
  end = clock();

  std::cout << "division took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

  bool isSame = true;

  for (long i = 0; i < SHE.getNumFactors(); i++) {
    isSame &= (decryption[i] == numerators[i]/denominator);
    if (decryption[i] != numerators[i]/denominator) {
      std::cout << "failed: " << numerators[i] << "/" << denominator << "!=" << decryption[i] << std::endl;
    } else {
      //std::cout << numerators[i] << "/" << denominator << "=" << decryption[i] << std::endl;
    }
  }

  return isSame;
}


bool testDivision() {
  std::cout << "oh boy.. :^)" << std::endl;
  long t = 257;
  NTL::ZZ q = NTL::GenPrime_ZZ(2348);
  //long d = 32768;
  //long d = 65536;
  long d = 4096;
  NTL::ZZ w = NTL::power2_ZZ(400);
  YASHE SHE(t,q,d,8,w);

  clock_t start, end;

  start = clock();

  NTL::ZZ_pX secretKey = SHE.keyGen();

  end = clock();

  std::cout << "Keygen took: " 
    << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

  srand(time(0));

  for (long i = 0; i < 100; i++) {

    long denominator, numerator;
    //do {
      denominator = rand() % t;
      numerator = rand() % t;
    //} while (denominator > numerator);

    YASHE_CT ciphertextN = SHE.encrypt(numerator);
    YASHE_CT ciphertextD = SHE.encrypt(denominator);

    start = clock();


    YASHE_CT::div(ciphertextN, ciphertextN, ciphertextD);

    end = clock();

    std::cout << "division took: " 
      << double(end - start)/CLOCKS_PER_SEC << "seconds" << std::endl;

    long output = SHE.decrypt(ciphertextN, secretKey);

    long desired;
    if (denominator != 0) {
      desired = numerator/denominator;
    } else {
      desired = numerator;
    }

    std::cout << "error: " << desired - output << "\t" << numerator << "/" << denominator << std::endl;

    //std::cout << numerator<<"/" << denominator << "= " << output << "(" << desired << ")" << std::endl;
  }

  return true;
}

void testReadWrite() {
  YASHE SHE = YASHE::readFromFile("8BitFHE");
  NTL::ZZ_pX sk = SHE.keyGen();
}


int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
