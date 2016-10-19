#ifndef FERRUM_TM_UTIL_HPP_
#define FERRUM_TM_UTIL_HPP_

#include "ferrum/minsky.hpp" // this must remain here
#include "ferrum/minsky.hpp"

#include <vector>

namespace ferrum {
  template <typename WCountType = int, typename RecordType = WCountType>
  inline RecordType
  populate_obs_lists(const minsky::CountList& doc_multi,
		     std::vector<int>& words_in_docs,
		     std::vector<RecordType>& word_type_counts);

}

#include "ferrum/tm_util.tcc"

#endif
