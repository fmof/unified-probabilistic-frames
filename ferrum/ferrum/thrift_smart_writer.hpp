#ifndef FERRUM_LIBNAR_THRIFT_SMART_WRITER_HPP_
#define FERRUM_LIBNAR_THRIFT_SMART_WRITER_HPP_

#include "ferrum/util.hpp"

namespace ferrum {
  namespace thrift {
    template <typename P>
    class ThriftSmartWriter : public ferrum::SmartWriter {
    public:
      ThriftSmartWriter<P>(const std::string& base);
      void reset_base(const std::string& base);
      virtual ~ThriftSmartWriter();
      virtual P* get();
      virtual P* get(const std::string& suffix) = 0;
      virtual ThriftSmartWriter<P>* clone() const = 0;

      virtual void advance(const std::string& suffix);
      template <typename Q, typename S> friend void save(ThriftSmartWriter<Q>* tsw, S& obj);
      virtual std::string form_id(const std::string& partial);
    protected:
      virtual void _save();
      boost::shared_ptr<P> protocol_;
    };

    template <typename P, typename S> void save(ThriftSmartWriter<P>* tsw, S& obj);
  }
}

#include "ferrum/thrift_smart_writer.tcc"

#endif
