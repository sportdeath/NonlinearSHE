#include <boost/lexical_cast.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <NTL/ZZ_p.h>

namespace boost {
  namespace serialization {
    template<class Archive>
      void save(Archive & ar, const NTL::ZZ_p & m, const unsigned int version) {
        std::string str(boost::lexical_cast<std::string>(rep(m)));
        ar & boost::serialization::make_nvp("str", str);
      }

    template<class Archive>
      void load(Archive & ar, NTL::ZZ_p & m, const unsigned int version) {
        std::string str;
        ar & boost::serialization::make_nvp("str", str);
        m = NTL::conv<NTL::ZZ_p>(NTL::to_ZZ(str.c_str()));
      }

    template<class Archive>
      void serialize(Archive & ar, NTL::ZZ_p & m, const unsigned int version) {
        boost::serialization::split_free(ar, m, version);
      }
  }
}

