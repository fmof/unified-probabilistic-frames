#ifndef FERRUM_ENTITIES_TCC_
#define FERRUM_ENTITIES_TCC_

#include "ferrum/words.hpp" // for VirtualVocabulary

#include <vector>

namespace ferrum {
  template <typename E>
  TemplatedEntityInterface<E>::TemplatedEntityInterface(const E& e) : entity(e) {
  }

  template <typename E>
  void TemplatedEntityInterface<E>::update_p_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const {
    vec->update_predicate_counts(vv, &entity, counts);
  }

  template <typename E>
  void TemplatedEntityInterface<E>::update_r_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const {
    vec->update_relation_counts(vv, &entity, counts);
  }
}

#endif
