#ifndef FERRUM_DATA_UTIL_TCC_
#define FERRUM_DATA_UTIL_TCC_

#include "ferrum/data_util.hpp"
#include "ferrum/timer.hpp"

namespace ferrum {
  namespace db {
    template <typename V, typename Conv>
    V Redis::list_item(const std::string& key, int idx, Conv converter) {
      RedisQuery rq(key);
      RedisQueryOptions opts;
      opts.lrange.start = idx;
      opts.lrange.end = idx;
      rq.options(opts);
      this->operator()(rq.lrange());
      const std::vector<std::string>& vals = rq.vec_value();
      if(vals.size() != 1) {
	ERROR << "Expected return of 1 element, got " << vals.size();
	V ret_val{};
	return ret_val;
      }
      V ret_val( converter(vals[0]) );
      return ret_val;
    }

    template <typename V, typename Conv>
    void Redis::list_item(const std::string& key, int idx, V* ret_val, Conv converter) {
      RedisQuery rq(key);
      RedisQueryOptions opts;
      opts.lrange.start = idx;
      opts.lrange.end = idx;
      rq.options(opts);
      this->operator()(rq.lrange());
      const std::vector<std::string>& vals = rq.vec_value();
      if(vals.size() != 1) {
	ERROR << "Expected return of 1 element, got " << vals.size();
      } else {
	converter(vals[0], ret_val);
      }
    }

    template <typename V>
    V PrimFromStr<V>::operator()(const std::string& v) {
      V ret{};
      std::stringstream ss;
      ss << v;
      ss >> ret;
      return ret;
    }

    template <typename V, typename Converter>
    V Redis::hget(const std::string& key, const std::string& field, Converter converter) {
      RedisQuery rq(key, field);
      rq.hget();
      this->operator()(rq);
      if(rq.is_set()) {
	V ret( (converter(rq.value())) );
	return ret;
      } else {
	return V{};
      }
    }

    // template <typename V, typename Converter>
    // int Redis::hset(const std::string& key, const std::string& field, const V& value) {
    //   std::string sval( (Converter(value)) );
    //   RedisQuery rq(key, field, sval);
    //   this->operator()(rq.hset());
    //   return std::stoi(rq.value());
    // }
  }
}

#endif
