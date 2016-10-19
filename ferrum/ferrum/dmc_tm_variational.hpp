#ifndef FERRUM_LIBNAR_DMC_TM_VARIATIONAL_HPP_
#define FERRUM_LIBNAR_DMC_TM_VARIATIONAL_HPP_

#include "ferrum/concrete.hpp"
#include "concrete_util/uuid_util.h"
#include "ferrum/minsky.hpp" // this must remain here
#include "ferrum/crtlda_variational.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/sage_defs.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/thrift_smart_writer.hpp"
#include "ferrum/util.hpp"

#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <thrift/protocol/TProtocol.h>
#include <time.h>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ferrum {
  /**
   * Variational inference for smoothed LDA.
   */
  class DMCTMVariational {
  private:
    typedef minsky::SimpleDoc DType;
    typedef UniformHyperSeedWeightedInitializer DiscreteVariationalInitializer;
    typedef std::vector< double > Vector1D;
    typedef Vector1D STopicType;
    typedef std::string VocabType;
    typedef std::vector<double> ModelTopicType;

    typedef std::vector< std::vector< double > > Vector2D;
    typedef std::pair< double, double > DoublePair;
    int num_docs_;

    Vector2D var_topic_usage_params_;
    // GLOBAL: Templates x Verbs
    std::vector< STopicType > var_topic_word_params_;

    // GLOBAL: Templates x Verbs
    Vector2D buffer_topic_word_params_;

    int num_topics_;
    int num_words_;

    Vector1D topic_usage_hypers_;
    Vector1D word_hypers_;

    // These are sparse views.
    std::vector< std::vector< int > > words_in_docs_;
    std::vector< std::vector< int > > word_type_counts_;

    DiscreteVariationalInitializer vi_init_;

    bool initialized_;

    template <typename CorpusT>
    CorpusT* corpus_for_labeling(BaseSituationLabeler*);

  public:
    typedef concrete::util::TCompactProtocol ThriftProtocol;
    DMCTMVariational();
    DMCTMVariational(int nt, const std::vector<double>& usage_h, const std::vector<double>& word_h);
    DMCTMVariational(int nt, int vsize, double usage_h, double word_h);
    ~DMCTMVariational();
    void num_topics(int nt);
    int num_topics();

    template <typename WCountType = int, typename RecordType = WCountType>
    RecordType populate_obs_lists(const minsky::CountList& doc_multi,
    				  std::vector<int>& words_in_docs,
    				  std::vector<RecordType>& word_type_counts);

    /**
     * Allocate and initialize global parameters
     */
    template <typename CorpusT>
    void init(const DiscreteVariationalInitializer&,
	      const Vocabulary<VocabType>&,
	      CorpusT* corpus // for polymorphism
	      );
    template <typename CorpusT>
    void init_batch(const CorpusT* corpus, // for polymorphism
		    const StringDiscreteKindPrinter& print_struct);

    //Vector2D get_distributions(const Vector2D* ptr);
    //Vector2D get_distributions(std::vector< STopicType >* vec);
    const Vector1D& buffer_topic_word_params(size_t i);
    void update_topic_assignments
    (
     Vector1D* assign_params,
     const Vector1D& usage_params,
     const Vector1D& expected_word_counts
     );
    void update_topic_usage(std::vector<double>* var_assign,
			    const std::vector<double>& buffer);
    void _update_topic_word_usage(int index_template, double i1 = 0.0, double i2 = 1.0);
    void update_topic_word_usage(int num_threads, double i1 = 0.0, double i2 = 1.0);
    
    // template <typename P>
    // void write_situations
    // (
    //  const DType& doc,
    //  ferrum::thrift::ThriftSmartWriter<P>* tsw,
    //  const Vector2D& t_assign,
    //  const Vector2D& s_assign
    //  );
    template <typename P>
    void write_usage_posterior
    (
     ferrum::thrift::ThriftSmartWriter<P>* tsw,
     const Vector1D& usage_post
     );

    template <typename CorpusT>
    double e_step
    (
     const VStrategy& strategy,
     int learn_iter,
     const size_t batch_start,
     const size_t batch_end,
     BaseSituationLabeler* sit_lab,
     const size_t offset = 0,
     int epoch = 0
     );

    void m_step(VStrategy& strategy, unsigned int batch_size);

    void update_hypers();

    void _print_in_learn
    (
     const VStrategy& strategy, 
     const StringDiscreteKindPrinter* const print_struct,
     DKVWriters* sw_wrapper,
     int epoch,
     int learn_iter,
     bool last_iter,
     bool model_changed
     );

    template <typename CorpusT, template <typename> class TSW, typename... TSWArgs>
    void learn(const CorpusT* corpus, // for polymorphism
	       VStrategy& strategy,
	       int epoch,
	       const StringDiscreteKindPrinter& print_struct,
	       DKVWriters* sw_wrapper,
	       TSWArgs... tsw_args);

    minsky::residual::ResidualTopicModel create_minsky_view(const ferrum::Vocabulary<std::string>& vocab);

  };
}

#include "ferrum/dmc_tm_variational.tcc"

#endif
