#ifndef FERRUM_DATA_THRIFT_PROTOCOL_DEFS_HPP_
#define FERRUM_DATA_THRIFT_PROTOCOL_DEFS_HPP_

#include <boost/shared_ptr.hpp>
#include <thrift/protocol/TProtocol.h>
#include <thrift/protocol/TBinaryProtocol.h>
#include <thrift/protocol/TCompactProtocol.h>
#include <thrift/protocol/TJSONProtocol.h>

namespace ferrum {
  namespace thrift {
    typedef apache::thrift::protocol::TBinaryProtocol TBinaryProtocol;
    typedef apache::thrift::protocol::TCompactProtocol TCompactProtocol;
    typedef apache::thrift::protocol::TJSONProtocol TJSONProtocol;
  }
}
#endif
