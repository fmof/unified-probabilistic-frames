#ifndef FERRUM_MINSKY_REDIS_CORPUS_HPP_
#define FERRUM_MINSKY_REDIS_CORPUS_HPP_

#include "ferrum/minsky.hpp"
#include "ferrum/redis_corpus.hpp"

#include <map>

namespace minsky {
  template <typename Doc>
  minsky::LabelMap<unsigned int>
  label_types(ferrum::RedisCorpus<Doc>* corpus);
}

#endif

#include "ferrum/minsky_redis_corpus.tcc"
