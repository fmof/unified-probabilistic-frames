#include "ferrum/crtlda_minsky.hpp"
#include "ferrum/logging.hpp"

namespace ferrum {
  const void* MinskyEntityCounter::inner_pred_vocab(const VirtualVocabulary* vv) const {
    return reinterpret_cast<const void* >(vv->downcast<std::string>() );
  }
  const void* MinskyEntityCounter::inner_rel_vocab(const VirtualVocabulary* vv) const {
    return reinterpret_cast<const void* >(vv->downcast<std::string>() );
  }
  MinskyEntityCounter::MinskyEntityCounter(const minsky::AnnotationLevel::type& al) : al_(al) {
  }
  MinskyEntityCounter::~MinskyEntityCounter() {
  }
  void MinskyEntityCounter::count_predicates(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const {
    for(const auto& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	if(pa.annot_level == al_) {
	  counts[pa.predicate.word] += 1;
	}
      }
    }
  }
  void MinskyEntityCounter::count_relations(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const {
    for(const auto& mention : entity.mentions) {
      for(const minsky::PredArg& pa : mention.structures) {
	if(pa.annot_level == al_) {
	  counts[pa.relation] += 1;
	}
      }
    }
  }
}
