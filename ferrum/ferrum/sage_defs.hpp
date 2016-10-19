#ifndef PUBLIC_UPF
#ifndef CRTLDA_SAGE_DEFS_H_
#define CRTLDA_SAGE_DEFS_H_

#include "ferrum/dmc.hpp"
#include "ferrum/lock.hpp"
#include "ferrum/loglin.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/minsky.hpp"
#include "omp.h"
#include "ferrum/optimize.hpp"
#include "ferrum/util.hpp"
#include "ferrum/wtm.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory> //shared_ptr
#include <ostream>
#include <stdlib.h>
#include <string>
#include <thread>
#include <time.h>
#include <unordered_set>
#include <utility>
#include <vector>


#ifndef SAGE_VERSION
#define SAGE_VERSION 1
#endif

namespace ferrum {
  enum SageTopicRegularization {
    L2 = 0,
    IMPROPER = 1
  };

  std::istream& operator>>(std::istream& in, SageTopicRegularization &how);
  double compute_sage_regularizer(double value, SageTopicRegularization how, double prev_val);
  double compute_grad_sage_regularizer(double value, SageTopicRegularization how, double prev_val);

  template <typename MaxentSupport, typename MaxentWeights>
  struct SageCClosure {
  public:
    typedef loglin::UnigramMaxent<MaxentSupport, MaxentWeights> MaxentModel;
    std::vector<double>* topic_counts;
    std::vector<double>* background;
    double marginal_count;
    MaxentModel* maxent;
    double regularizer_multiplier;
    SageTopicRegularization regularizer_type = SageTopicRegularization::L2;
    std::vector<MaxentSupport>* sparse_which_items = NULL;
    std::vector<int>* sparse_which_counts = NULL;
    std::vector<double>* prev_eta = NULL;
    bool is_model_ready = false;
  };

  struct SageCClosureIntVector : public SageCClosure<int, std::vector<double> > {
  };

  enum TopicInitializerChoice {
    UNIFORM = 0,
    SUBSET = 1,
    BACKGROUND = 2
  };

  void complete_background_init(int num_words, std::vector<double>* raw_bck, double min_log_prob = -200, double min_freq = 0.0001);

  template <typename EtaType>
  class SageTopic {
  public:
    typedef std::vector<double> TauType;
    typedef int MaxentSupportType;
    typedef loglin::UnigramMaxent<MaxentSupportType, EtaType> MaxentModel;
    typedef SageCClosureIntVector ClosureType;
  protected:
    std::shared_ptr< std::vector<double> > background_;
    std::unique_ptr< EtaType > eta_;
    //    TauType tau_;
    std::vector<double>* summed_weights_;
    
    // This model shares the eta_ weights
    MaxentModel maxent_;

    double tau_hyper_;
    double log_partition_;
    int support_size_;

    SageTopicRegularization regularization_type_;
    optimize::OptimizationMethod opt_how_;
    std::shared_ptr< optimize::Minimizer > minimizer_;

    double sparsity_thres_;
    bool sparse_;
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
    void __eta();
  public:
    typedef EtaType Eta;

    SageTopic<EtaType>();
    SageTopic<EtaType>
    (
     int support_size,
     double kappa, 
     SageTopicRegularization reg_type
     );
    SageTopic<EtaType>
    (
     int support_size,
     double kappa, 
     SageTopicRegularization reg_type,
     optimize::OptimizationMethod opt_meth
     );
    SageTopic<EtaType>(const SageTopic<EtaType>& other);
    SageTopic<EtaType>(SageTopic<EtaType>&& other);
    SageTopic<EtaType>& operator=(const SageTopic<EtaType>& other);
    SageTopic<EtaType>& operator=(SageTopic<EtaType>&& other);
    virtual ~SageTopic();
    double& operator[](size_t idx);
    void prepare();
    void push_to_maxent();
    MaxentModel* maxent_model();
    std::vector<double>* background();
    void background(const std::shared_ptr<std::vector<double> >& background);
    double log_partition();
    //void renormalize(const std::vector<double>& logterm_sums);
    void renormalize();
    void sparse(bool b);
    bool sparse();
    void sparsity_threshold(double d);
    double sparsity_threshold();
    minsky::Frame create_minsky_frame(int which_voc);
    /**
     * Copy the parameter values.
     * Recompute the un-normalized, log-weights,
     * "store" them in the underlying Maxent model,
     * and renormalize.
     */
    void eta(const EtaType& eta);
    /**
     * Return a copy of the stored (residual) parameters
     */
    EtaType eta();
    double l_probability(const MaxentSupportType& obj);
    double probability(const MaxentSupportType& obj);

    template <typename OutputType> OutputType as();

    template <typename OutputType> OutputType eta_as(bool normalize = true);
    template <typename OutputType> OutputType eta_as(bool normalize = true) const;

    static double elbo_static(ClosureType *fparams, const double *weights, const int n, bool renorm);
    double elbo(const std::vector<double>& counts, bool renorm);
    double elbo(const std::vector<int>& counts, const std::vector<MaxentSupportType>& items, bool renorm);

    // these next two satisfy the general optimize interface
    static double elbo_eval(void *fparams, const double* point, const int size);
    static void elbo_grad(void *fparams, const double* point, double* grad, const int size);

    // this satisfies the lbfgs interface
    static double liblbfgs_elbo(void *fparams, const lbfgsfloatval_t *trial_weights,
				lbfgsfloatval_t *lbfgs_grad, const int n, const lbfgsfloatval_t step);
    std::shared_ptr<optimize::Function> get_optimizable_func(ClosureType* params);
    int fit_topic(const std::vector<double>* counts, double reg_mult);
    int fit_topic_svi(const std::vector<double>* counts, double reg_mult, double i1, double i2);
    double variance_first_deriv(const std::pair<double,double>& pair, double e_i);
    double variance_second_deriv(const std::pair<double,double>& pair, double e_i);
    void update_variances(std::vector<std::pair<double,double> >* params, double i1, double i2);

    // This gets an initial point.
    // Following the original SAGE implementation, this initializes to the zero vector
    // (which makes sense, because in expectation, due to the sparsity-inducing prior, 
    // it should be zero).
    Eta get_optimization_initial_point();
    int support_size();

    minsky::Frame create_minsky_frame();
  };

  std::istream& operator>>(std::istream& in, TopicInitializerChoice &how);

  struct SageInitializer {
  private:
    int num_words_;
    int num_threads_;
    TopicInitializerChoice fit_topic_how_;
    SageTopicRegularization reg_topic_how_;
    int num_docs_for_init_;
    struct {
      double threshold_;
      bool opt_;
      bool serialize_;
    } sparsity_info;
    struct {
      unsigned int topic_subset;
      std::shared_ptr<std::mutex> subset_mut;
    } num_calls;

    /**
     * Provide parameters eta over V elements for a SageTopic maxent model. 
     * This effectively returns "-background + signed random noise", 
     * in order to produce a SAGE topic 
     *   p(w | eta) \propto exp(background + eta)
     * that is near uniform.
     */
    std::vector<double> topic_uniform(const std::vector<double>& hyper,
				      std::shared_ptr<std::vector<double> > background,
				      double num_types);
    /**
     * Provide parameters eta over V elements for a SageTopic maxent model. 
     * This produces signed random noise around 0,
     * in order to produce a SAGE topic 
     *   p(w | eta) \propto exp(background + eta)
     * that is close to the background distribution. 
     * Effectively, this results in a maximum-likelihood topic 
     * **distribution**.
     */
    std::vector<double> topic_background(const std::vector<double>& hyper,
					 std::shared_ptr<std::vector<double> > background,
					 double num_types);
    /**
     * Provide parameters eta over V elements for a SageTopic maxent model. 
     * This randomly selects d' documents in the corpus (or all D if 
     * d' > D), where 
     *   d' = 10000/(average number of tokens per doc).
     * This fits a SAGE topic 
     *   p(w | eta) \propto exp(background + eta)
     * where the background parameters are fixed and shared across
     * all topics.
     */
    template <typename CorpusT>
    std::vector<double> topic_fit_subset(const std::vector<double>& hyper,
					 CorpusT* corpus,
					 std::shared_ptr< std::vector<double> > background,
					 minsky::MinskyAnnotationWrapper maw
					 //minsky::AnnotationLevel::type annot_level,
					 //minsky::StructureType::type st_level
					 );

    /**
     * Provide parameters eta over V elements for a SageTopic maxent model. 
     */
    template <typename CorpusT, typename Vocab>
    std::vector<double> topic(const std::vector<double>& hyper,
			      CorpusT* corpus,
			      const Vocab* vocab,
			      std::shared_ptr<std::vector<double> > background,
			      minsky::MinskyAnnotationWrapper maw
			      );
    std::vector<std::pair<double, double> > gamma_variance(size_t size);
  public:
    SageInitializer();
    SageInitializer(int num_words);
    SageInitializer(int num_words, int num_threads);
    SageInitializer(int num_words, int num_thread,
		    TopicInitializerChoice fit_how, 
		    SageTopicRegularization reg_topic_how);
    
    const SageTopicRegularization& regularize_how() const;
    /**
     * Return the maximum number of threads that *may* be
     * used during initialization. Not all initializations are
     * multithreaded, so the number of threads actually used 
     * could be significantly less than the value returned here.
     */
    int num_threads(bool safe = false);
    /**
     * TopicInitializerChoice::SUBSET, set the number of documents
     * to draw. If not set, then the following number are drawn:
     *
     *            10000 * # documents
     *          ----------------------
     *            # tokens in corpus
     *
     */
    void num_docs_for_init(int);

    void sparse(bool thrift, bool opt, double thresh);

    /**
     * Provide multinomial parameters over K elements, were each
     * p_k \propto 1/K + U(-1/K, 1/K)
     */
    std::vector<double> assignment(size_t num);

    double threshold();
    bool sparse();

    template <typename CorpusT, typename Vocab>
    std::vector<double> words(const std::vector<double>& hyper,
			      CorpusT* corpus,
			      const Vocab* vocab,
			      std::shared_ptr<std::vector<double> > background) {
      minsky::MinskyAnnotationWrapper maw;
      maw.word_level = minsky::WordAnnotation::ORTHOGRAPHIC;
      return topic(hyper, corpus, vocab, background, maw);
    }

    template <typename CorpusT, typename Vocab>
    std::vector<double> verbs(const std::vector<double>& hyper,
			      CorpusT* corpus,
			      const Vocab* vocab,
			      std::shared_ptr<std::vector<double> > background) {
      minsky::MinskyAnnotationWrapper maw;
      maw.annot_level = minsky::AnnotationLevel::SYNTAX;
      maw.st_level = minsky::StructureType::PREDICATE;
      return topic(hyper, corpus, vocab, background, maw);
    }

    template <typename CorpusT, typename Vocab>
    std::vector<double> arcs(const std::vector<double>& hyper,
			     CorpusT* corpus,
			     const Vocab* vocab,
			     std::shared_ptr<std::vector<double> > background) {
      minsky::MinskyAnnotationWrapper maw;
      maw.annot_level = minsky::AnnotationLevel::SYNTAX;
      maw.st_level = minsky::StructureType::RELATION;
      return topic(hyper, corpus, vocab, background, maw);
    }

    /**
     * Provide Dirichlet parameters gamma over K elements.
     * Given N words in a **document** and initial parameters alpha, 
     * each gamma_k is
     *   gamma_k = alpha_k + M/K + U(-M/K, M/K)
     * This is a better option for large corpora, vs
     * usage(const std::vector<double>&).
     */
    std::vector<double> usage(const std::vector<double>& hyper, const double N);
    std::vector<double> usage_templates(const std::vector<double>& hyper, const double N);

    /**
     * Provide Dirichlet parameters gamma over K elements.
     * Given M words in a corpus, and initial parameters alpha, 
     * each gamma_k is
     *   gamma_k = alpha_k + M/K + U(-M/K, M/K)
     * For large corpora, this results in parameters that yield
     * very flat Dirichlet draws. In such cases, it is better to use
     * usage(const std::vector<double>&, const int).
     */
    std::vector<double> usage(const std::vector<double>& hyper);

    std::vector<double> slots(const std::vector<double>& hyper);

    std::vector<std::pair<double, double> > verb_variance(size_t size);
    std::vector<std::pair<double, double> > arc_variance(size_t size);
  };

  struct SageStrategy {
    SageStrategy() {
    }
    double em_frobenius_threshold = 1E-6;
    double eta_density_threshold = 1E-4;
    int num_learn_iters = 100;
    int num_e_iters = 25;
    int num_m_iters = 1;
    int num_e_threads = 1;
    int num_m_threads = 1;
    int hyper_update_min = 20;
    int hyper_update_iter = 5;
    int update_model_every = 5;
    int partial_restarts = 0;
    int num_learn_restart_iters = 25;
    int num_e_restart_iters = 25;
    int print_topics_every = 5;
    int print_topics_k = 10;
    int print_usage_every = 5;
    int em_verbosity = 1;
    bool heldout = false;
    bool never_update_model = false;
  };
}

#include "ferrum/sage_defs.tcc"

#endif
#endif
