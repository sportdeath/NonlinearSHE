#include <boost/lexical_cast.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/binary_object.hpp>

#include <NTL/ZZ.h>

namespace boost {
  namespace serialization {
    template<class Archive>
      void save(Archive & ar, const NTL::ZZ & m, const unsigned int version) {
        std::string str(boost::lexical_cast<std::string>(m));
        //unsigned char * bytes;
        //std::size_t size = NumBytes(m);
        //bytes = (unsigned char *) malloc (size);
        //std::cout << "size: " << size << std::endl;
        //std::cout << m << std::endl;

        //BytesFromZZ(bytes, m, size);
        ar & boost::serialization::make_nvp("str", str);
            //boost::serialization::make_binary_object(bytes, size));
        //std::cout << "mmmm" << std::endl;
      }

    template<class Archive>
      void load(Archive & ar, NTL::ZZ & m, const unsigned int version) {
        //unsigned char * bytes;
        //std::size_t size;

        //std::cout << "ohh" << std::endl;

        //ar & boost::serialization::make_nvp("data", 
            //boost::serialization::make_binary_object(bytes, size));
        //std::cout << "size: " << size << std::endl;
        //ZZFromBytes(m, bytes, size);
        std::string str;
        ar & boost::serialization::make_nvp("str", str);
        m = NTL::to_ZZ(str.c_str());
      }

    template<class Archive>
      void serialize(Archive & ar, NTL::ZZ & m, const unsigned int version) {
        boost::serialization::split_free(ar, m, version);
      }
  }
}

