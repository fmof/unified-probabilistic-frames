#ifndef __CRTLDA_MINSKY_THRIFT_H_
#define __CRTLDA_MINSKY_THRIFT_H_

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/minsky.hpp"

namespace ferrum {
  // an entity counter for minsky::Entity
  class MinskyEntityCounter : public VirtualEntityCounterCRTP< minsky::Entity > {
  protected:
    const void* inner_pred_vocab(const VirtualVocabulary* vv) const;
    const void* inner_rel_vocab(const VirtualVocabulary* vv) const;
  public:
    typedef minsky::Entity E;
    minsky::AnnotationLevel::type al_;
    MinskyEntityCounter(const minsky::AnnotationLevel::type& al);
    virtual ~MinskyEntityCounter();
    void count_predicates(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const;
    void count_relations(const VirtualVocabulary* vv, const E& entity, std::vector<double>& counts) const;
  };

  template <typename CCorpus, typename BC>
  class BDocBackgroundCounter {
  public:
    BDocBackgroundCounter();
    BDocBackgroundCounter& in_memory_as(std::shared_ptr<std::vector<double> >);
    BDocBackgroundCounter& defined_by(const VirtualVocabulary*);
    //BDocBackgroundCounter& how(const MinskyEntityCounter*);
    BDocBackgroundCounter& with(BC*);
    BDocBackgroundCounter& over(CCorpus*);

    template <typename VocabElemType>
    void compute_background(BackgroundComputeParams bcp = BackgroundComputeParams());

  protected:
    std::shared_ptr<std::vector<double> > counts_;
    BC* bc_;
    const VirtualVocabulary* vv_;
    //const MinskyEntityCounter* mec_;
    CCorpus* corpus_;
  }; //BDocBackgroundComputer

  template <typename CCorpus, typename BC>
  class EDocBackgroundCounter {
  public:
    EDocBackgroundCounter();
    EDocBackgroundCounter& in_memory_as(std::shared_ptr<std::vector<double> >);
    EDocBackgroundCounter& defined_by(const VirtualVocabulary*);
    EDocBackgroundCounter& how(const MinskyEntityCounter*);
    EDocBackgroundCounter& with(BC*);
    EDocBackgroundCounter& over(CCorpus*);

    template <typename VocabElemType>
    void compute_predicate_background(BackgroundComputeParams bcp = BackgroundComputeParams());

    template <typename VocabElemType>
    void compute_relation_background(BackgroundComputeParams bcp = BackgroundComputeParams());
  protected:
    std::shared_ptr<std::vector<double> > counts_;
    BC* bc_;
    const VirtualVocabulary* vv_;
    const MinskyEntityCounter* mec_;
    CCorpus* corpus_;
  }; //EDocBackgroundComputer

}

#endif

#include "ferrum/crtlda_minsky.tcc"
