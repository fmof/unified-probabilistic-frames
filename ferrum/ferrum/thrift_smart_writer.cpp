#include "ferrum/thrift_protocol_defs.hpp"
#include "ferrum/thrift_smart_writer.hpp"

namespace ferrum {
  namespace thrift {
    template <typename P>
    ThriftSmartWriter<P>::ThriftSmartWriter(const std::string& base) : 
      ferrum::SmartWriter(base),
      protocol_(NULL) {
    }

    template <typename P>
    ThriftSmartWriter<P>::~ThriftSmartWriter() {
    }

    template <typename P>
    void ThriftSmartWriter<P>::reset_base(const std::string& base) {
      this->base_ = base;
    }

    template <typename P>
    std::string ThriftSmartWriter<P>::form_id(const std::string& partial) {
      return partial;
    }

    template <typename P>
    P* ThriftSmartWriter<P>::get() {
      return this->get("");
    }

    template <typename P> void ThriftSmartWriter<P>::_save() {
    }
  }
}

template class ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TBinaryProtocol>;
template class ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TCompactProtocol>;
template class ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TJSONProtocol>;
