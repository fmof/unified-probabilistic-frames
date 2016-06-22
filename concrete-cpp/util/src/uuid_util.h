#ifndef CONCRETECPP_UTIL_UUID_UTIL_H_
#define CONCRETECPP_UTIL_UUID_UTIL_H_

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include <boost/random/mersenne_twister.hpp>

namespace concrete { 
  namespace util {
    class uuid_factory {
    private:
      boost::uuids::basic_random_generator<boost::mt19937> gen;
    public:
      uuid_factory() {
      };
      std::string get_uuid() {
	boost::uuids::uuid uuid = gen();
	return boost::uuids::to_string(uuid);
      };
    };
  }
}

#endif
