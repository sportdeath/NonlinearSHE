#include <vector>
#include <ctime>
#include <cmath>

#include <NTL/ZZ.h>

#include <YASHE/YASHE.hpp>
#include <YASHE/cipherText.hpp>
#include <YASHE/functions.hpp>

#include "gtest/gtest.h"

class YASHE8BitTest : public ::testing::Test {
  protected:

    static long d, t, sigma;
    static NTL::ZZ q, w;
    static NTL::ZZ_pX secretKey;
    static YASHE SHE;

    virtual void SetUp() {
      srand(time(0));
      NTL::ZZ_p::init(q);
    }

};

long YASHE8BitTest::t = 257;
NTL::ZZ YASHE8BitTest::q = NTL::GenPrime_ZZ(800);
long YASHE8BitTest::d = 512; // This parameter is set low so tests are faster
                             // it must be ~ 20000 for a security parameter of 128
                             // 22016 = 2^9 * 43 works well as it has 5376
                             // distinct factors - giving it a batch of that size.
long YASHE8BitTest::sigma = 8;
NTL::ZZ YASHE8BitTest::w = NTL::power2_ZZ(200);
YASHE YASHE8BitTest::SHE = YASHE(t,q,d,sigma,w);
NTL::ZZ_pX YASHE8BitTest::secretKey = SHE.keyGen();



TEST_F(YASHE8BitTest, RandomKeyPolyBounds) {
  NTL::ZZ_pX randPoly = SHE.randomKeyPoly();
  for (long i = 0; i <= deg(randPoly); i++ ) {
    ASSERT_TRUE(randPoly[i] == 0 |
                randPoly[i] == 1 |
                randPoly[i] == -1);
  }
}

TEST_F(YASHE8BitTest, RandomKeyPolyDifferent) {
  NTL::ZZ_pX randPoly1 = SHE.randomKeyPoly();
  NTL::ZZ_pX randPoly2 = SHE.randomKeyPoly();

  ASSERT_NE(randPoly1, randPoly2);
}

TEST_F(YASHE8BitTest, RandomErrPolyBounds) {
  NTL::ZZ_pX randPoly = SHE.randomKeyPoly();

  for (long i = 0; i <= deg(randPoly); i++ ) {
    ASSERT_TRUE(rep(randPoly[i]) < sigma * 6 |
                rep(randPoly[i]) > q - 6 * sigma);
  }
}


TEST_F(YASHE8BitTest, RandomERRPolyDifferent) {
  NTL::ZZ_pX randPoly1 = SHE.randomErrPoly();
  NTL::ZZ_pX randPoly2 = SHE.randomErrPoly();

  ASSERT_NE(randPoly1, randPoly2);
}


TEST(YASHETest, RadixDecompInt) {
  YASHE SHE(2,NTL::ZZ(257),5,8,NTL::ZZ(3));

  std::vector<NTL::ZZ> decomp;
  SHE.radixDecomp(decomp, NTL::ZZ(200));
  
  std::vector<long> result {2, 0, 1, 1, 2, 0};

  for (long i = 0; i < decomp.size(); i++) {
    ASSERT_EQ(decomp[i], result[i]);
  }
}


TEST(YASHETest, RadixDempPoly) {
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

  for (long i = 0; i < 6; i++) {
    for (long j = 0; j < 5; j++) {
      ASSERT_EQ(decomp[i][j], outs[i][j]);
    }
  }
}

TEST(YASHETest, PowersOfRadix) {
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

  for (long i = 0; i < 6; i++) {
    for (long j = 0; j < 5; j++) {
      ASSERT_EQ(powers[i][j], outs[i][j]);
    }
  }
}

TEST_F(YASHE8BitTest, DecompPowers) {
  NTL::ZZ_pX testPoly1 = SHE.randomErrPoly();
  NTL::ZZ_pX testPoly2 = SHE.randomErrPoly();

  NTL::ZZ_pXModulus cycloMod(NTL::conv<NTL::ZZ_pX>(SHE.getCycloModX()));

  NTL::ZZ_pX mulPoly;
  NTL::MulMod(mulPoly, testPoly1, testPoly2, cycloMod);

  std::vector<NTL::ZZ_pX> powers, decomp;
  SHE.powersOfRadix(powers, testPoly1);
  SHE.radixDecomp(decomp, testPoly2);

  NTL::ZZ_pX dotPoly;
  SHE.dot(dotPoly, powers, decomp);

  ASSERT_EQ(mulPoly, dotPoly);
}

TEST_F(YASHE8BitTest, EncryptDecrypt) {
  long message = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  for (long i = 0; i <=SHE.getMaxDegree(); i++) {
    ASSERT_EQ(message, decryption);
  }
}


TEST_F(YASHE8BitTest, AdditionByConstant) {
  long message = rand() % t;
  long constant = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);

  YASHE_CT::add(ciphertext, ciphertext, constant);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  ASSERT_EQ((message + constant) % t, decryption);
}


TEST_F(YASHE8BitTest, MultiplyByConstant) {
  long message = rand() % t;
  long constant = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);

  YASHE_CT::mul(ciphertext, ciphertext, constant);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  ASSERT_EQ((message * constant) % t, decryption);
}


TEST_F(YASHE8BitTest, Add) {
  long message1 = rand() % t;
  long message2 = rand() % t;

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);
  long decryption = SHE.decrypt(ciphertext1, secretKey);

  ASSERT_EQ((message1 + message2) % t, decryption);
}


TEST_F(YASHE8BitTest, Subtract) {
  long message1 = rand() % t;
  long message2 = rand() % t;

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  YASHE_CT::sub(ciphertext1, ciphertext1, ciphertext2);
  long decryption = SHE.decrypt(ciphertext1, secretKey);

  ASSERT_EQ((message1 + t - message2) % t, decryption);
}

TEST_F(YASHE8BitTest, Multiply) {
  long message1 = rand() % t;
  long message2 = rand() % t;

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);

  long decryption = SHE.decrypt(ciphertext1, secretKey);

  ASSERT_EQ((message1 * message2) % t, decryption);
}


TEST_F(YASHE8BitTest, EvaluatePolynomial) {
  long degree = t - 1;
  NTL::ZZ_pX poly;
  poly.SetLength(degree + 1);
  for (long i = 0; i <= degree; i++) {
    poly[i] = rand() % t;
  }

  poly[degree] = 1;

  long message = rand() % t;

  YASHE_CT ciphertext = SHE.encrypt(message);
  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  long result = 0;
  for (long i = degree - 1; i >= 0; i--) {
    result = (result * message + rep(poly[i])) % t;
  }

  ASSERT_EQ(result, decryption);
}


TEST_F(YASHE8BitTest, DivisionByConstant) {
  long message = rand() % t;
  long constant = rand() % t;

  NTL::ZZ_pX poly = 
    Functions::functionToPoly(Functions::divideByConstant(constant), t);

  YASHE_CT ciphertext = SHE.encrypt(message);
  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);

  long decryption = SHE.decrypt(ciphertext, secretKey);

  ASSERT_EQ(message/constant, decryption);
}


TEST_F(YASHE8BitTest, Division) {
  long message1 = rand() % t;
  long message2 = rand() % t;

  YASHE_CT ciphertext1 = SHE.encrypt(message1);
  YASHE_CT ciphertext2 = SHE.encrypt(message2);

  YASHE_CT::div(ciphertext1, ciphertext1, ciphertext2);

  long decryption = SHE.decrypt(ciphertext1, secretKey);

  ASSERT_NEAR(message1/message2, decryption, 5);
}


TEST_F(YASHE8BitTest, EncryptDecryptBatch) {
  std::vector<long> message(SHE.getNumFactors());

  for (long i = 0; i <= SHE.getNumFactors(); i++) {
    message[i] = rand() % t;
  }

  YASHE_CT ciphertext = SHE.encryptBatch(message);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++) {
    ASSERT_EQ(message[i], decryption[i]);
  }
}


TEST_F(YASHE8BitTest, AddBatch) {
  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = rand() % t;
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++)
    ASSERT_EQ((message1[i] + message2[i]) % t, decryption[i]);
}

TEST_F(YASHE8BitTest, SubtractBatch) {
  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = rand() % t;
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  YASHE_CT::sub(ciphertext1, ciphertext1, ciphertext2);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++)
    ASSERT_EQ((message1[i] + t - message2[i]) % t, decryption[i]);
}

TEST_F(YASHE8BitTest, MultiplyBatch) {
  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());
  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = rand() % t;
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++)
    ASSERT_EQ((message1[i] * message2[i]) % t, decryption[i]);
}


TEST_F(YASHE8BitTest, DivisionByConstantBatch) {
  std::vector<long> message(SHE.getNumFactors());

  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message[i] = rand() % t;
  }

  YASHE_CT ciphertext = SHE.encryptBatch(message);

  long constant = 0;
  while (constant == 0) {
   constant = rand() % t;
  }

  NTL::ZZ_pX poly = 
    Functions::functionToPoly(Functions::divideByConstant(constant), t);

  YASHE_CT::evalPoly(ciphertext, ciphertext, poly);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++)
    ASSERT_EQ(message[i]/constant , decryption[i]);
}


TEST_F(YASHE8BitTest, DivisionBatch) {
  std::vector<long> message1(SHE.getNumFactors());
  std::vector<long> message2(SHE.getNumFactors());

  for (long i = 0; i < SHE.getNumFactors(); i++) {
    message1[i] = rand() % t;
    message2[i] = 0;
    while (message2[i] == 0) {
      message2[i] = rand() % t;
    }
  }

  YASHE_CT ciphertext1 = SHE.encryptBatch(message1);
  YASHE_CT ciphertext2 = SHE.encryptBatch(message2);

  YASHE_CT::div(ciphertext1, ciphertext1, ciphertext2);

  std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

  ASSERT_EQ(decryption.size(), SHE.getNumFactors());

  for (long i = 0; i < decryption.size(); i++)
    ASSERT_NEAR(message1[i]/message2[i] , decryption[i], 5);
}

int main(int argc, char ** argv) {
  //::testing::GTEST_FLAG(filter) = "*DecompPowers*";
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
