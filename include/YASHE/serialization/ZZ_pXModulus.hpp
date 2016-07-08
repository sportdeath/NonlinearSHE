#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <NTL/ZZ_pX.h>

namespace boost {
  namespace serialization {
    template<class Archive>
      void serialize(Archive & ar, NTL::ZZ_pXModulus & m, const unsigned int version) {
        ar & m.f;
        ar & m.UseFFT;
        ar & m.n;
        ar & m.k;
        ar & m.FRep;
        ar & m.HRep;
      }
  }
}

