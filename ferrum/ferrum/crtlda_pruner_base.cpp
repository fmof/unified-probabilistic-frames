#include "ferrum/crtlda_pruner_base.hpp"

namespace ferrum {
  GeneralizedDependencyPath::GeneralizedDependencyPath() :
    arc(),
    dependency_index(-1) {
  }
  GeneralizedDependencyPath::~GeneralizedDependencyPath() {
  }

  PathGeneralizer::PathGeneralizer
  (
   const std::vector<concrete::TaggedToken>* pos_list,
   const std::vector<concrete::Dependency>* dep_list,
   unsigned int max_hops) : 
    agg_list_(),
    plp_(pos_list),
    dlp_(dep_list),
    max_hops_(max_hops) {
  }
  bool PathGeneralizer::generalize_to_v(size_t curr_dep_idx, unsigned int curr_hop) {
    bool able = false;
    if(curr_hop > max_hops_) return able;
    const concrete::Dependency& curr_dep = dlp_->operator[](curr_dep_idx);
    if(curr_dep.gov < 0) return able;
    seen_dep_idx_.insert(curr_dep_idx);
    GeneralizedDependencyPath path;
    switch(plp_->operator[](curr_dep.gov).tag[0]) {
    case 'V': // verb
      path.arc = curr_dep.edgeType;
      path.dependency_index = curr_dep_idx;
      agg_list_.push_back(path);
      able = true;
      break;
    case 'J': // adjective
      {
	const size_t dsize = dlp_->size();
	for(size_t i = 0; i < dsize; ++i) {
	  const concrete::Dependency& next_dep = dlp_->operator[](i);
	  if(next_dep.dep != curr_dep.gov) continue;
	  if(i == curr_dep_idx) continue;
	  if(seen_dep_idx_.count(i)) continue;
	  // note: do NOT add it here
	  bool child_able = generalize_to_v(i, curr_hop+1);
	  able |= child_able;
	  if(child_able) break;
	}
	if(able) {
	  path.arc = curr_dep.edgeType;
	  path.dependency_index = curr_dep_idx;
	  agg_list_.push_back(path);
	}
      }
      break;
    default: // nothing
      break;
    }
    return able;
  }

  const std::list<GeneralizedDependencyPath>& PathGeneralizer::get_path() {
    return agg_list_;
  }

  std::string PathGeneralizer::form_gov_tagging(const std::vector<concrete::TaggedToken>& ttl) {
    std::stringstream ss;
    size_t idx = 0;
    for(const GeneralizedDependencyPath& gdp : agg_list_) {
      int gov_tok_idx = dlp_->operator[](gdp.dependency_index).gov;
      ss << ttl[gov_tok_idx].tag;
      ++idx;
      if(idx < agg_list_.size()) {
	ss << ":";
      }
    }
    return ss.str();
  }

  std::string PathGeneralizer::form_dep_tagging(const std::vector<concrete::TaggedToken>& ttl) {
    std::stringstream ss;
    size_t idx = 0;
    for(const GeneralizedDependencyPath& gdp : agg_list_) {
      int dep_tok_idx = dlp_->operator[](gdp.dependency_index).dep;
      ss << ttl[dep_tok_idx].tag;
      ++idx;
      if(idx < agg_list_.size()) {
	ss << ":";
      }
    }
    return ss.str();
  }

  std::string PathGeneralizer::form_arc_tagging() {
    std::stringstream ss;
    size_t idx = 0;
    for(const GeneralizedDependencyPath& gdp : agg_list_) {
      ss << dlp_->operator[](gdp.dependency_index).edgeType;
      ++idx;
      if(idx < agg_list_.size()) {
	ss << ":";
      }
    }
    return ss.str();
  }

  size_t PathGeneralizer::size() {
    return agg_list_.size();
  }
}
