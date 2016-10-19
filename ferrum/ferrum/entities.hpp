#ifndef FERRUM_ENTITIES_HPP_
#define FERRUM_ENTITIES_HPP_

#include "ferrum/words.hpp" // for VirtualVocabulary

#include <vector>

namespace ferrum {
  // forward declare
  template <typename E>
  class TemplatedEntityInterface;
  
  // forward declare
  class VirtualEntityCounter;

  class EntityInterface {
  public:
    template <typename E> const TemplatedEntityInterface<E>* downcast() const {
      //return dynamic_cast< const TemplatedEntityInterface<E>* >(this);
      return reinterpret_cast< const TemplatedEntityInterface<E>* >(this);
    }
    virtual void update_p_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const = 0;
    virtual void update_r_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const = 0;
  };

  template <typename E>
  class TemplatedEntityInterface : public EntityInterface {
  public:
    E entity;
    TemplatedEntityInterface<E>(const E& e);
    void update_p_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const;
    void update_r_count(const VirtualEntityCounter* vec, const VirtualVocabulary* vv, std::vector<double>& counts) const;
  };
  
  class VirtualEntityCounter {
  protected:
    virtual const void* inner_pred_vocab(const VirtualVocabulary* vv) const = 0;
    virtual const void* inner_rel_vocab(const VirtualVocabulary* vv) const = 0;
  public:
    virtual ~VirtualEntityCounter() {
    }
    virtual void update_predicate_counts(const VirtualVocabulary* vv, const void *ei, std::vector<double>& counts) const = 0;
    virtual void update_relation_counts(const VirtualVocabulary* vv, const void *ei, std::vector<double>& counts) const = 0;

    template <typename W>
    unsigned int serialize_string(const Vocabulary<W>* vocab, const std::vector<double>& counts, std::stringstream& ss) const {
      //      const Vocabulary< W >* vocab = vv->downcast<W>();
      unsigned int expected_size = 0;
      for(size_t i = 0; i < counts.size(); ++i) {
	const W& word = vocab->word(i);
	std::string countstring = std::to_string( counts[i] );
	ss << word << " " << countstring << " ";
	expected_size += word.size() + countstring.size() + 2;
      }
      return expected_size;
    }

    template <typename E>
    TemplatedEntityInterface<E> make_entity_interface(const E& ent) const {
      return TemplatedEntityInterface<E>(ent);
    }

    template <typename W>
    const Vocabulary<W>* predicate_vocab(const VirtualVocabulary* vv) const {
      return reinterpret_cast<const Vocabulary<W>*>( inner_pred_vocab(vv) );
    }
    template <typename W>
    const Vocabulary<W>* relation_vocab(const VirtualVocabulary* vv) const {
      return reinterpret_cast<const Vocabulary<W>*>( inner_rel_vocab(vv) );
    }
  };

  template <typename E>
  class VirtualEntityCounterCRTP : public VirtualEntityCounter {
  public:
    virtual ~VirtualEntityCounterCRTP() {
    }
    virtual void count_predicates(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const = 0;
    virtual void count_relations(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const = 0;

    void update_predicate_counts(const VirtualVocabulary* vv, const void* entityp, std::vector<double>& counts) const {
      const E& entity = *(reinterpret_cast<const E*>(entityp));
      count_predicates(vv, entity, counts);
    }
    void update_relation_counts(const VirtualVocabulary* vv, const void* entityp, std::vector<double>& counts) const {
      const E& entity = *(reinterpret_cast<const E*>(entityp));
      count_relations(vv, entity, counts);
    }
  };
}

namespace ferrum {
  template <typename E>
  using TemplatedEntityInterface = ferrum::TemplatedEntityInterface<E>;

  typedef ferrum::VirtualEntityCounter VirtualEntityCounter;

  template <typename E>
  using VirtualEntityCounterCRTP = ferrum::VirtualEntityCounterCRTP<E>;
}

#endif

#include "ferrum/entities.tcc"
