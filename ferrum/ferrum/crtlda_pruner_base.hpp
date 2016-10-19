#ifndef FERRUM_LIBNAR_CRTLDA_PRUNER_BASE_H_
#define FERRUM_LIBNAR_CRTLDA_PRUNER_BASE_H_

#include "ferrum/dmc.hpp"
#include "ferrum/util.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/words.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include <time.h>

// for pair
#include "map"
#include <utility>
#include <unordered_set>
#include <string>
#include <vector>

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

namespace ferrum {
  struct Toolnames {
    std::string dep_parse_tool = "col-ccproc-deps";
    std::string entity_mention_tool = "Stanford";
    std::string situation_mention_tool = "Semafor";
  };

  class GeneralizedDependencyPath {
  public:
    std::string arc;
    int dependency_index;
    GeneralizedDependencyPath();
    ~GeneralizedDependencyPath();
  };

  class PathGeneralizer {
  public:
    PathGeneralizer(const std::vector<concrete::TaggedToken>* pos_list, const std::vector<concrete::Dependency>* dep_list, unsigned int max_hops);
    bool generalize_to_v(size_t curr_dep_idx, unsigned int curr_hop);
    std::string form_gov_tagging(const std::vector<concrete::TaggedToken>& ttl);
    std::string form_dep_tagging(const std::vector<concrete::TaggedToken>& ttl);
    std::string form_arc_tagging();
    size_t size();
    const std::list<GeneralizedDependencyPath>& get_path();
  private:
    std::list<GeneralizedDependencyPath> agg_list_;
    const std::vector<concrete::TaggedToken>* plp_;
    const std::vector<concrete::Dependency>* dlp_;
    unsigned int max_hops_;
    std::set<size_t> seen_dep_idx_;
  };

  template <typename M>
  class EntityMentionPruner {
  protected:
    const concrete::Communication& communication;
  public:
    EntityMentionPruner<M>(const concrete::Communication& comm) : communication(comm) {
    }
    virtual std::vector<M> prune(const concrete::EntityMention& conc_mention) const = 0;
  };
}
#endif
