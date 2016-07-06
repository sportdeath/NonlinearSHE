#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

namespace boost {
  namespace serialization {
    template<class Archive>
      void save(Archive & ar, const NTL::ZZ_pX & m, const unsigned int version) {
        // copy to std vector of zz
        std::vector<NTL::ZZ_p> coeffiecients(deg(m) + 1);
        for (long i = 0; i <= deg(m); i++) {
          coeffiecients[i] = m[i];
        }
        ar & boost::serialization::make_nvp("vector", coeffiecients);
      }

    template<class Archive>
      void load(Archive & ar, NTL::ZZ_pX & m, const unsigned int version) {
        std::vector<NTL::ZZ_p> coeffiecients;
        ar & boost::serialization::make_nvp("vector", coeffiecients);
        m.SetLength(coeffiecients.size());
        for (long i = 0; i < coeffiecients.size(); i++) {
          m[i] = coeffiecients[i];
        }
      }

    template<class Archive>
      void serialize(Archive & ar, NTL::ZZ_pX & m, const unsigned int version) {
        boost::serialization::split_free(ar, m, version);
      }
  }
}

