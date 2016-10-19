#ifndef FERRUM_MINSKY_REDIS_CORPUS_TCC_
#define FERRUM_MINSKY_REDIS_CORPUS_TCC_

#include "ferrum/minsky.hpp"
#include "ferrum/minsky_redis_corpus.hpp"
#include "ferrum/redis_corpus.hpp"

#include <map>

namespace minsky {
  template <typename Doc>
  minsky::LabelMap<unsigned int>
  label_types(ferrum::RedisCorpus<Doc>* corpus) {
    minsky::LabelMap<unsigned int> lmap;
    for(typename ferrum::RedisCorpus<Doc>::const_iterator it = corpus->begin();
	it != corpus->end();
	++it) {
      const Doc& doc = *(it->document);
      for(const minsky::Label& label : doc.labels) {
	lmap[label] += 1;
      }
    }
    return lmap;
  }
}

#endif
