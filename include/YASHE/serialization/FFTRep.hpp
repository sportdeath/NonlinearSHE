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
      void serialize(Archive & ar, NTL::FFTRep & m, const unsigned int version) {
        ar & m.k;
        ar & m.MaxK;
        //ar & make_array<long> (m.tbl, 4);
        ar & m.NumPrimes;
      }
  }
}

