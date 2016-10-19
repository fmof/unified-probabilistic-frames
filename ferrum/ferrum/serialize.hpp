#ifndef FERRUM_SERIALIZE_H_
#define FERRUM_SERIALIZE_H_

#include "ferrum/logging.hpp"

// for serialization
#include <fstream>
#include <iostream>

// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace serialize {
template < typename Serializable >
inline void gzip_serialize(Serializable* object, std::string to_file, bool on_condition = true) {
  if(on_condition) {
    std::ofstream ofs(to_file, std::ios::out|std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::gzip_compressor());
    out.push(ofs);
    boost::archive::binary_oarchive oa(out);
    oa << (*object);
    INFO << "see " << to_file << " for serialized object.";
  }
}

}
#endif
