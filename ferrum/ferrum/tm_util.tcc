#ifndef FERRUM_TM_UTIL_TCC_
#define FERRUM_TM_UTIL_TCC_

#include "ferrum/minsky.hpp"

#include <vector>

namespace ferrum {
  template <typename WCountType, typename RecordType>
  inline
  RecordType
  populate_obs_lists(const minsky::CountList& doc_multi,
		     std::vector<int>& wid,
		     std::vector<RecordType>& wtc) {
    RecordType num_words = (RecordType)0;
    for(const auto& count_pair : doc_multi.icounts) {
      const int word = count_pair.first;
      wid.push_back(word);
      RecordType c = (RecordType)(count_pair.second);
      wtc.push_back(c);
      num_words += c;
    }
    return num_words;
  }
}

#endif
