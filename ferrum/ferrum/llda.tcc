#ifndef FERRUM_LIBNAR_LLDA_TCC_
#define FERRUM_LIBNAR_LDA_TCC_

#include "ferrum/llda.hpp"
#include "ferrum/tm_util.hpp"
#include "ferrum/timer.hpp"
#include <utility>

namespace ferrum {
  template <typename CorpusT>
  CorpusT* DMCTMVariational::corpus_for_labeling(BaseSituationLabeler* sit_lab) {
    CorpusT* corp_for_labeling = NULL;
    if(sit_lab != NULL && sit_lab->perform_labeling()) {
      corp_for_labeling = sit_lab->corpus< CorpusT* >();
      assert(corp_for_labeling != NULL);
      //doc = corp_for_labeling->operator[](di + offset);
    }
    return corp_for_labeling;
  }
  template <typename CorpusT>
  void DMCTMVariational::init
  (
   const DiscreteVariationalInitializer& vi,
   const Vocabulary<VocabType>& vocab,
   CorpusT* corpus // for polymorphism
   ) {
    ferrum::Timer init_time(__func__);
    this->vi_init_ = vi;
    num_words_ = vocab.num_words();
    var_topic_word_params_.resize(num_topics_);
    buffer_topic_word_params_.resize(num_topics_);
    for(int t = 0; t < num_topics_; ++t) {
      //STopicType v_topic(num_words_, 0.0);
      STopicType v_topic{
	vi_init_.words(word_hypers_)
      };
      var_topic_word_params_[t] = std::move(v_topic);
      buffer_topic_word_params_[t] = Vector1D(num_words_, 0.0);
    }
    initialized_ = true;
    INFO << "Initialization complete";
  }

  template <typename CorpusT>
  void DMCTMVariational::init_batch
  (
   const CorpusT* corpus,
   const StringDiscreteKindPrinter& print_struct
   ) {
    num_docs_ = corpus->num_docs();
    // ensure that all document-specific globals are empty
    var_topic_usage_params_.clear();
    words_in_docs_.clear();
    word_type_counts_.clear();
    // and now initialize
    int num_words = 0;
    for(typename CorpusT::const_iterator doc_it = corpus->begin();
	doc_it != corpus->end();
	++doc_it) {
      const DType& doc = *(doc_it->document);
      size_t di = doc_it->iteration_idx;
      minsky::CountList bow = minsky::bow_counts(doc);
      const size_t num_words_in_doc = minsky::num_words(bow);
      var_topic_usage_params_.push_back(vi_init_.usage_template(topic_usage_hypers_, num_words_in_doc));
      words_in_docs_.push_back(std::vector<int>());
      word_type_counts_.push_back(std::vector<int>());
      num_words += populate_obs_lists(bow, words_in_docs_[di], word_type_counts_[di]);
    }
  }

  template <typename CorpusT, template <typename> class TSW, typename... TSWArgs >
  void DMCTMVariational::learn
  (
   const CorpusT* corpus,
   VStrategy& strategy,
   int epoch,
   const StringDiscreteKindPrinter& print_struct,
   DKVWriters* sw_wrapper,
   TSWArgs... tsw_args
   ) {
    bool last_iter = false;
    bool force_update = false;
    //bool model_changed = false;
    if(strategy.em_verbosity >= 0) {
      INFO << "EM Epoch " << epoch;
    }
    const size_t total_num_docs = corpus->num_docs();
    const size_t batch_size =
      strategy.batch_size <= 0 ? total_num_docs : (size_t)strategy.batch_size;
    const size_t num_batches = 
      (size_t)((total_num_docs % batch_size) ?
	       (total_num_docs / batch_size + 1) :
	       (total_num_docs / batch_size));
    size_t corpus_start = 0, 
      corpus_end = batch_size;
    INFO << "Running inference with " << num_batches << " batches.";
    // break the provided subset into different batches
    for(size_t batch_i = 0; batch_i < num_batches; ++batch_i) {
      ferrum::Timer batch_timer(std::string(__func__) + " :: batch " + std::to_string(batch_i));
      // Initialize the provided subset
      INFO << "Batch " << batch_i << ": preparing corpus of from documents [" << corpus_start << ", " << corpus_end << ")";
      CorpusT subset_corpus = corpus->subset(corpus_start, corpus_end);
      init_batch(&subset_corpus, print_struct);
      size_t num_loaded = num_docs_;
      //bool label_docs = (strategy.label_every > 0);
      // ConcreteSituationLabeler<ThriftProtocol, TSW> sit_lab(label_docs);
      // sit_lab.template make_with_args< std::string, TSWArgs... >("situations.tcompact", tsw_args...);
      // if(label_docs) {
      // 	sit_lab.corpus(&subset_corpus);
      // }

      double batch_elbo =
	e_step<CorpusT>
	(
	 strategy,
	 epoch,
	 0, num_loaded, // the start/end coordinates are a bit useless (right now), since we're saying that each batch is just a subset of the *entire* (passed-in) corpus. We could instead define batches *within* subsets of the passed-in corpus. In that case, any of the following could change: the start point, the end point, or the offset (below)
	 NULL, // &sit_lab,
	 0, // the batch/subset offset
	 epoch
	 );
      if(strategy.compute_elbo || strategy.elbo_as_break) {
	INFO << "Epoch " << epoch << ", from [" << corpus_start << ", "<< corpus_end << "), e-step elbo is " << batch_elbo;
      }
      if(! strategy.heldout) {
	m_step(strategy, (unsigned int)batch_size);
      }	
      if( strategy.hyper_update_iter > 0 ) {
	update_hypers();
      }
      // if(! strategy.never_update_model) {
      // 	if(corpus_end % strategy.update_model_every == 0 || force_update || last_iter) {
      // 	  model_changed = true;
      // 	}
      // }
      //_print_in_learn(strategy, &print_struct, sw_wrapper, epoch, corpus_end, last_iter, model_changed);
      corpus_start = corpus_end + 1;
      corpus_end += batch_size;
      if(corpus_end >= total_num_docs) {
	corpus_end = total_num_docs;
	last_iter = true;
      }
    }
  }
}
#endif
