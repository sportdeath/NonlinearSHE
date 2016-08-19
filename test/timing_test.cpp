#include <vector>
#include <ctime>
#include <cmath>

#include <NTL/ZZ.h>

#include <YASHE/YASHE.hpp>
#include <YASHE/cipherText.hpp>
#include <YASHE/functions.hpp>

#include "gtest/gtest.h"

class YASHETiming : public ::testing::TestWithParam<std::vector<long>> {
  public:

    const static long numTrials = 2;

    NTL::ZZ_pX secretKey;
    YASHE SHE;

    clock_t start, end;
    double totalTime;

    virtual void SetUp() {
      long t = GetParam()[0];
      long log2Q = GetParam()[1];
      long d = GetParam()[2];
      long log2W = GetParam()[3];

      srand(time(0));

      long sigma = 8;

      //NTL::ZZ q = NTL::GenPrime_ZZ(log2Q);
      //NTL::ZZ w = NTL::power2_ZZ(log2W);

      //NTL::ZZ_p::init(q);

      SHE = YASHE(t, log2Q, d, sigma, log2W);

      secretKey = SHE.keyGen();

      resetTimer();
    }

    virtual void TearDown() {
      printTimingResults();
    }

    YASHE_CT getRandomCiphertext(
        std::vector<long> & messages, 
        bool includeZero = true) {

      messages.resize(SHE.getNumFactors());

      for (long i = 0; i < SHE.getNumFactors(); i++) {
        if (includeZero) {
          messages[i] = rand() % SHE.getPModulus();
        } else {
          do {
            messages[i] = rand() % SHE.getPModulus();
          } while (messages[i] == 0);
        }
      }

      return SHE.encryptBatch(messages);
    }

    double clocksToMS() {
      return double(end - start)/(CLOCKS_PER_SEC/1000.);
    }

    void printTimingResults() {
      double timePerBatch = totalTime/numTrials;
      double timePerOperation = timePerBatch/SHE.getNumFactors();
      std::cout << "Total Time: " << totalTime << " ms"
                << "\t Time/Batch: " << timePerBatch  << " ms"
                << "\t Time/Operation: " << timePerOperation << " ms"
                << std::endl;
    }

    void resetTimer() {
      totalTime = 0;
    }

    void startTimer() {
      start = clock();
    }

    void pauseTimer() {
      end = clock();
      totalTime += clocksToMS();
    }

};

TEST_P(YASHETiming, MultiplicativeDepth) {
  std::cout << "log2q: " << NTL::NumBits(SHE.getCModulus()) << std::endl;

  //for (long trial = 0; trial < numTrials; trial ++) {
    
    std::vector<long> messages1, messages2;
    YASHE_CT ciphertext1 = getRandomCiphertext(messages1);
    YASHE_CT ciphertext2 = getRandomCiphertext(messages2);

    long level = 0;
    while(true) {
      std::cout << "At level " << level << std::endl;
      YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);
      std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

      bool isCorrect = true;

      for (long i = 0; i < SHE.getNumFactors(); i++) {
        messages1[i] = (messages1[i] * messages2[i]) % SHE.getPModulus();
        isCorrect &= 
          (messages1[i] == decryption[i]);
      }

      if (isCorrect) {
        level += 1;
      } else {
        break;
      }
    }

    std::cout << "Got to level: " << level << std::endl;
  //}
}

//TEST_P(YASHETiming, Encryption) {
  //for (long trial = 0; trial < numTrials; trial ++) {
    //std::vector<long> messages;
    //startTimer();
    //YASHE_CT ciphertext = getRandomCiphertext(messages);
    //pauseTimer();
  //}
//}

//TEST_P(YASHETiming, Decryption) {
  //for (long trial = 0; trial < numTrials; trial ++) {
    //std::vector<long> messages;
    //YASHE_CT ciphertext = getRandomCiphertext(messages);

    //startTimer();
    //std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);
    //pauseTimer();

    //for (long i = 0; i < SHE.getNumFactors(); i++)
      //ASSERT_EQ(messages[i], decryption[i]);
  //}
//}

//TEST_P(YASHETiming, Addition) {
  //for (long trial = 0; trial < numTrials; trial ++) {
    //std::vector<long> messages1, messages2;
    //YASHE_CT ciphertext1 = getRandomCiphertext(messages1);
    //YASHE_CT ciphertext2 = getRandomCiphertext(messages2);

    //startTimer();
    //YASHE_CT::add(ciphertext1, ciphertext1, ciphertext2);
    //pauseTimer();

    //std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

    //for (long i = 0; i < SHE.getNumFactors(); i++)
      //ASSERT_EQ((messages1[i] + messages2[i]) % SHE.getPModulus(), decryption[i]);
  //}
//}

//TEST_P(YASHETiming, Multiplication) {
  //for (long trial = 0; trial < numTrials; trial ++) {
    //std::vector<long> messages1, messages2;
    //YASHE_CT ciphertext1 = getRandomCiphertext(messages1);
    //YASHE_CT ciphertext2 = getRandomCiphertext(messages2);

    //startTimer();
    //YASHE_CT::mul(ciphertext1, ciphertext1, ciphertext2);
    //pauseTimer();

    //std::vector<long> decryption = SHE.decryptBatch(ciphertext1, secretKey);

    //for (long i = 0; i < SHE.getNumFactors(); i++)
      //ASSERT_EQ((messages1[i] * messages2[i]) % SHE.getPModulus(), decryption[i]);
  //}
//}

TEST_P(YASHETiming, RandomPoly) {
  for (long trial = 0; trial < numTrials; trial ++) {
    std::vector<long> messages;
    YASHE_CT ciphertext = getRandomCiphertext(messages);

    long polyDegree = SHE.getPModulus() - 1;
    std::vector<long> poly(polyDegree);
    for (long deg = 0; deg < polyDegree; deg++) {
      poly[deg] = rand() % SHE.getPModulus();
    }

    startTimer();
    YASHE_CT::evalPoly(ciphertext, ciphertext, poly);
    pauseTimer();

    std::vector<long> decryption = SHE.decryptBatch(ciphertext, secretKey);

    for (long i = 0; i < SHE.getNumFactors(); i++) {
      long result = 0;
      for (long deg = polyDegree - 1; deg >= 0; deg--) {
        result = (result * messages[i] + poly[deg]) % SHE.getPModulus();
      }
      ASSERT_EQ(result, decryption[i]);
    }
  }
}


std::vector< std::vector<long> > values8Bit {
  {   257,  583, 22016,  32},
  {   257,  437, 22016,  32},
  {   257,  389, 22016,  32},
  {   257, 1195, 66048,  32},
  {   257, 1055, 66048,  32},
  {   257,  791, 66048,  32},
};

std::vector< std::vector<long> > values16Bit {
  { 65537,  992, 32768,  32},
  { 65537,  878, 32768,  32},
  { 65537,  658, 32768,  32},
  { 65537, 1850, 66536,  32},
  { 65537, 1627, 66536,  32},
  { 65537, 1218, 66536,  32}
};

std::vector< std::vector<long> > valuesOptimized {
  //{   257,  400, 22016,  32},
  //{   257,  584, 22016,  192},
  //{   257,  500, 22016,  192},
  {   257,  1195, 66048,  512},
};

INSTANTIATE_TEST_CASE_P(
    OptimizedValues,
    YASHETiming, 
    ::testing::ValuesIn(valuesOptimized)
    );

//INSTANTIATE_TEST_CASE_P(
    //8BitTests,
    //YASHETiming, 
    //::testing::ValuesIn(values8Bit)
    //);

//INSTANTIATE_TEST_CASE_P(
    //16BitTests,
    //YASHETiming, 
    //::testing::ValuesIn(values16Bit)
    //);

 //division
 //equality

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
