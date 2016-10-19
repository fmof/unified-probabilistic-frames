#ifndef FERRUM_LIBNAR_THRIFT_SMART_WRITER_TCC_
#define FERRUM_LIBNAR_THRIFT_SMART_WRITER_TCC_

#include "ferrum/thrift_smart_writer.hpp"
#include "ferrum/unused_var.hpp"

template <typename P, typename S>
void ferrum::thrift::save(ThriftSmartWriter<P>* tsw, S& obj) {
  obj.write( tsw->protocol_.get() );
  tsw->_save();
}

namespace ferrum {
  namespace thrift {
    template <typename P>
    void ThriftSmartWriter<P>::advance(const std::string& suffix) {
      ignore_result( this->get(suffix));
    }
  }
}

#endif
