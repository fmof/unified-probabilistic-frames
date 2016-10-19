#ifndef __MINSKY_THRIFT_H_
#define __MINSKY_THRIFT_H_

#include "minsky/bdoc_types.h"
#include "minsky/edoc_types.h"
#include "minsky/frames_types.h"
#include "minsky/residual_global_slots_types.h"
#include "minsky/residual_tm_types.h"
#include "minsky/vocab_types.h"

#include <map>

namespace minsky {
  struct MinskyAnnotationWrapper {
    minsky::AnnotationLevel::type annot_level = minsky::AnnotationLevel::UNSPECIFIED;
    minsky::StructureType::type st_level    = minsky::StructureType::UNSPECIFIED;
    minsky::WordAnnotation::type  word_level  = minsky::WordAnnotation::UNSPECIFIED;
  };
  bool operator<(const MinskyAnnotationWrapper& x, const MinskyAnnotationWrapper& y);
  size_t num_entities(const EDoc& doc);
  std::map<int, int> predicate_histogram(const Entity& entity, size_t which_structure = 0);
  std::map<int, int> relation_histogram(const Entity& entity, size_t which_structure = 0);
  std::map<int, int> predicate_histogram_on_level(const Entity& entity, const AnnotationLevel::type& al);
  std::map<int, int> relation_histogram_on_level(const Entity& entity, const AnnotationLevel::type& al);
  std::map<int, int> predicate_histogram(const EDoc& doc, const AnnotationLevel::type& al);
  std::map<int, int> relation_histogram(const EDoc& doc, const AnnotationLevel::type& al);
  int find_level(const Mention& ent, const AnnotationLevel::type& al, bool print = true);
  minsky::Frame frame_unnormalized(const std::vector<double>& v);

  minsky::CountList& operator+=(minsky::CountList& recept, const minsky::CountList& other);
  minsky::CountList bow_counts(const SimpleDoc& doc);
  size_t num_words(const CountList& counts);
  size_t num_words(const SimpleDoc& doc);
  bool empty(const WordsClause& wc);

  template <typename Derived, typename DocType>
  class MinskyDocVocabUpdater {
  public:
    MinskyDocVocabUpdater<Derived, DocType>() {
    }
    virtual ~MinskyDocVocabUpdater() {
    }
    void operator()(DocType& doc, const std::vector<std::vector<int> >& vocab_mappers) {
      static_cast<Derived*>(this)->operator()(doc, vocab_mappers);
    }
  };
  class SyntacticEDocVocabUpdater : public MinskyDocVocabUpdater<SyntacticEDocVocabUpdater, minsky::EDoc> {
  public:
    SyntacticEDocVocabUpdater();
    void operator()(minsky::EDoc& doc, const std::vector<std::vector<int> >& vocab_mappers);
  };
  class PredArgEDocVocabUpdater : public MinskyDocVocabUpdater<PredArgEDocVocabUpdater, minsky::EDoc> {
  public:
    /**
     * The order of the vocab levels (per group, predicate then relation)
     */
    PredArgEDocVocabUpdater(const std::vector<minsky::MinskyAnnotationWrapper>& orders);
    void operator()(minsky::EDoc& doc, const std::vector<std::vector<int> >& vocab_mappers);
  private:
    //std::vector<minsky::MinskyAnnotationWrapper> orders_;
    std::map<minsky::MinskyAnnotationWrapper, int> orders_;
  };

  class SimpleDocVocabUpdater : public MinskyDocVocabUpdater<SimpleDocVocabUpdater, minsky::SimpleDoc> {
  public:
    SimpleDocVocabUpdater();
    void operator()(minsky::SimpleDoc& doc, const std::vector<std::vector<int> >& vocab_mappers);
  };

  struct LabelCmp {
    bool operator() (const minsky::Label& lhs, const minsky::Label& rhs) const;
  };

  template <typename Number>
  using LabelMap = std::map<minsky::Label, Number, LabelCmp>;

#ifdef THRIFT091
  std::ostream& operator<<(std::ostream&, const minsky::Label&);
#endif
  std::string get_label(const minsky::SimpleDoc&);
}

namespace ferrum {
  std::map<int, int> doc_multinomial(const minsky::EDoc& doc,
				     const minsky::MinskyAnnotationWrapper& maw);
  std::map<int, int> doc_multinomial(const minsky::SimpleDoc& doc,
				     const minsky::MinskyAnnotationWrapper& maw);
				     //const minsky::AnnotationLevel::type& annot_level,
				     //const minsky::StructureType::type& st_level);
  const std::vector<minsky::Entity>& get_entities(const minsky::EDoc& doc);
  const std::vector<minsky::Mention>& get_mentions(const minsky::Entity& entity);
  int num_entities(const minsky::EDoc& doc);
  int num_entities(const minsky::Entity& entity);
  int get_num_mentions(const minsky::EDoc& doc);
}

#endif
