#ifndef __CRTLDA_MINSKY_THRIFT_TCC_
#define __CRTLDA_MINSKY_THRIFT_TCC_

namespace ferrum {
  template <typename CCorpus, typename BC>
  BDocBackgroundCounter<CCorpus, BC>::BDocBackgroundCounter() {
  }

  template <typename CCorpus, typename BC>
  BDocBackgroundCounter<CCorpus, BC>& BDocBackgroundCounter<CCorpus, BC>::in_memory_as(std::shared_ptr<std::vector<double> > c) {
    counts_ = c;
    return *this;
  }
  template <typename CCorpus, typename BC>
  BDocBackgroundCounter<CCorpus, BC>& BDocBackgroundCounter<CCorpus, BC>::defined_by(const VirtualVocabulary* vv) {
    vv_ = vv;
    return *this;
  }
  template <typename CCorpus, typename BC>
  BDocBackgroundCounter<CCorpus, BC>& BDocBackgroundCounter<CCorpus, BC>::with(BC* bc) {
    bc_ = bc;
    return *this;
  }
  template <typename CCorpus, typename BC>
  BDocBackgroundCounter<CCorpus, BC>& BDocBackgroundCounter<CCorpus, BC>::over(CCorpus* corp) {
    corpus_ = corp;
    return *this;
  }

  template <typename CCorpus, typename BC>
  template <typename VocabElemType>
  void BDocBackgroundCounter<CCorpus, BC>::compute_background
  (
   BackgroundComputeParams bcp
   ) {
    bc_->before();
    typedef typename CCorpus::const_iterator const_iterator;
    typedef typename CCorpus::DocType DocType;
    for(const_iterator it = corpus_->begin(); it != corpus_->end(); ++it) {
      const DocType& document = *(it->document);
      minsky::CountList cl = minsky::bow_counts(document);
      for(const std::pair<int, int> pcount : cl.icounts) {
	counts_->operator[](pcount.first) += pcount.second;
      }
    }
    bc_->complete_background_compute(counts_.get(), bcp);
    const Vocabulary<VocabElemType>* wvocab = vv_->downcast<VocabElemType>();
    bc_->after(wvocab, corpus_, counts_.get(), "words");
  }

  ///////////////////////////////////////////////////////////////////////

  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>::EDocBackgroundCounter() {
  }

  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>& EDocBackgroundCounter<CCorpus, BC>::in_memory_as(std::shared_ptr<std::vector<double> > c) {
    counts_ = c;
    return *this;
  }
  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>& EDocBackgroundCounter<CCorpus, BC>::defined_by(const VirtualVocabulary* vv) {
    vv_ = vv;
    return *this;
  }
  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>& EDocBackgroundCounter<CCorpus, BC>::how(const MinskyEntityCounter* mec) {
    mec_ = mec;
    return *this;
  }
  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>& EDocBackgroundCounter<CCorpus, BC>::with(BC* bc) {
    bc_ = bc;
    return *this;
  }
  template <typename CCorpus, typename BC>
  EDocBackgroundCounter<CCorpus, BC>& EDocBackgroundCounter<CCorpus, BC>::over(CCorpus* corp) {
    corpus_ = corp;
    return *this;
  }

  template <typename CCorpus, typename BC>
  template <typename VocabElemType>
  void EDocBackgroundCounter<CCorpus, BC>::compute_predicate_background
  (
   BackgroundComputeParams bcp
   ) {
    bc_->before();
    typedef typename CCorpus::const_iterator const_iterator;
    typedef typename CCorpus::DocType DocType;
    for(const_iterator it = corpus_->begin(); it != corpus_->end(); ++it) {
      const DocType& document = *(it->document);
      for(const auto& entity : get_entities(document) ) {
	TemplatedEntityInterface< decltype(entity) > tei = mec_->make_entity_interface< decltype(entity) >(entity);
	tei.update_p_count(mec_, vv_, *counts_);
      }
    }
    bc_->complete_background_compute(counts_.get(), bcp);
    const Vocabulary<VocabElemType>* wvocab = mec_->predicate_vocab<VocabElemType>(vv_);
    bc_->after(wvocab, corpus_, counts_.get(), "predicate");
  }

  template <typename CCorpus, typename BC>
  template <typename VocabElemType>
  void EDocBackgroundCounter<CCorpus, BC>::compute_relation_background(BackgroundComputeParams bcp) {
    typedef typename CCorpus::const_iterator const_iterator;
    typedef typename CCorpus::DocType DocType;
    bc_->before();
    for(const_iterator it = corpus_->begin(); it != corpus_->end(); ++it) {
      const DocType& document = *(it->document);
      for(const auto& entity : get_entities(document) ) {
	TemplatedEntityInterface< decltype(entity) > tei = mec_->make_entity_interface< decltype(entity) >(entity);
	//es->update_relation_counts(vocab, &tei, counts);
	tei.update_r_count(mec_, vv_, *counts_);
      }
    }
    bc_->complete_background_compute(counts_.get(), bcp);
    const Vocabulary<VocabElemType>* wvocab = mec_->relation_vocab<VocabElemType>(vv_);
    bc_->after(wvocab, corpus_, counts_.get(), "relation");
  }


}
#endif
