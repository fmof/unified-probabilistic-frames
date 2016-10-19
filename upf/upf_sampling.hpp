#ifndef UPF_SAMPLING_H_
#define UPF_SAMPLING_H_

#include "ferrum/crtlda_defs.hpp"
#include "ferrum/crtlda_sampling.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/util.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_writers.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <time.h>

// for pair
#include <map>
#include <utility>
#include <unordered_set>
#include <string>
#include <vector>

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

namespace ferrum {
  /**
   * A collapsed Gibbs sampler for a discrete model.
   * This assumes that each observation is generated by a
   * categorical (discrete) distribution that is endowed
   * with a Dirichlet prior.
   */
  class CollapsedGibbsDMC {
  public:
    typedef ferrum::Vocabulary<std::string> VocabT;
    typedef minsky::EDoc DocT;
    typedef concrete::util::TCompactProtocol ThriftProtocol;
    CollapsedGibbsDMC(int num_templates, int num_slots, const ferrum::SymmetricHyperparams& shp);
    ~CollapsedGibbsDMC();
    void use_lexical(bool b);
    bool use_lexical();
    template <typename CorpusT> void init(int num_templates, 
					  const std::vector<int>& num_slots,
					  CorpusT* corpus,
					  const VocabT& gov_vocab,
					  const VocabT& rel_vocab);
    template <typename CorpusT> void init(int num_templates, int num_slots,
					  CorpusT* corpus,
					  const VocabT& gov_vocab,
					  const VocabT& rel_vocab);
    template <typename CorpusT> void init(CorpusT* corpus,
					  const VocabT& gov_vocab,
					  const VocabT& rel_vocab);
    // template <typename G, typename R> void heldout_init(int num_templates, 
    // 							const std::vector<int>& num_slots,
    // 							const VocabT& gov_vocab,
    // 							const VocabT& rel_vocab);
    //void heldout_init();
    
    int sample_template(const int curr_slot_val);
    int sample_slot(const int curr_template_val);
    void adjust_for_slot(int templ_val, int slot_val, int di,
			 int ei, int by);
    void adjust_for_template(int val, int di, int ei, int by);
    void unassign_template(int prev, int di, int ei);
    void assign_template(int sampled, int di, int ei);
    void unassign_slot(int templ, int prev, int di, int ei);
    void assign_slot(int templ, int sampled, int di, int ei);
    void update_floating_hists(const minsky::Entity& entity);
    std::vector< std::vector<double> > doc_template_params();
    void sampling_strategy(SamplingStrategy* ss);
    void print_dist_diag(const ferrum::StringDiscreteKindPrinter& print_struct, size_t K);
    template <typename CorpusT> void print_labeling(CorpusT* corpus, const std::string& f_name, const VocabT& gv, const VocabT& rv);

    template <typename CorpusT, template <typename> class TSW, typename... TSWArgs >
    void learn
    (
     CorpusT* corpus,
     int epoch,
     const StringDiscreteKindPrinter& print_struct,
     DKVWriters* sw_wrapper,
     bool heldout,
     TSWArgs... tsw_args
     );

    minsky::residual::ResidualUniqueSlots create_minsky_view(const VocabT& gov_vocab, const VocabT& rel_vocab) const;
    
    template <typename Protocol>
    void write_usage_posterior
    (
     ferrum::thrift::ThriftSmartWriter<Protocol>* tsw,
     const std::vector<double>& usage
     );
  private:
    int num_docs_;
    SamplingStrategy* sample_strategy_;
    SymmetricHyperparams shp_;

    int num_templates_;
    std::vector<int> num_slots_;

    // assignment arrays
    std::vector< std::vector<int> > template_assignments_;
    std::vector< std::vector<int> > slot_assignments_;
    typedef std::vector<int>* floating_counter_ptr;

    // count tables
    // num_docs * num_templates: how many times a given template is used in each document
    std::vector<std::vector<int> > c_template_;
    floating_counter_ptr c_doc_template_ptr_;
    // num_templates * num_govs: how many times a given template is used for a governor
    std::vector<std::vector<int> > c_template_gov_;
    // num_templates: how many times a template occurs
    std::vector<int>  c_template_sum_;
    // num_templates * num_slots_per_template: how many times a slot is used (per template) and by doc
    std::vector< std::vector<int> > c_slot_;
    // num_templates * num_slots_per_template * num rels
    std::vector<std::vector< std::vector<int> > > c_slot_rel_;

    gsl_rng_type *which_gsl_rng;
    gsl_rng *rnd_gen_;
    bool own_rng_;

    std::map<int, int> floating_gov_hist_;
    std::map<int, int> floating_rel_hist_;

    dmc::gdmc t_usage_dmc_;
    dmc::mfgdmc rel_dmc_;
    dmc::gdmc gov_dmc_;
    dmc::gdmc slot_dmc_;

    bool use_lexical_;
    bool use_lexical_set_;

    minsky::AnnotationLevel::type annot_level_ = minsky::AnnotationLevel::SYNTAX;
  };

  ////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////

  /**
   * A collapsed Gibbs sampler for a discrete model.
   * This assumes that each observation is generated by a
   * categorical (discrete) distribution that is endowed
   * with a Dirichlet prior.
   */
  class CollapsedGibbsDMCWithKinds {
  public:
    // syntactic frame-level vocab
    typedef ferrum::Vocabulary<std::string> VocabT;
    // semantic frame-level vocab
    typedef ferrum::Vocabulary<std::string> VocabF;
    typedef minsky::EDoc DocT;
    typedef concrete::util::TCompactProtocol ThriftProtocol;
    void latent_kinds(bool b);
    void use_lexical(bool b);
    bool use_lexical();
    CollapsedGibbsDMCWithKinds(int nt, int ns, int nf, int nr,
			       const SymmetricHyperparams& shp,
			       bool latent_kinds);
    ~CollapsedGibbsDMCWithKinds();

    void rnd_gen(gsl_rng* rng);

    template <typename CorpusT>
    void init(int num_templates, 
	      const std::vector<int>& num_slots,
	      CorpusT* corpus,
	      const VocabF& gov_kind_vocab,
	      const VocabF& rel_kind_vocab,
	      const VocabT& gov_vocab,
	      const VocabT& rel_vocab);
    template <typename CorpusT>
    void init(CorpusT* corpus,
	      const VocabF& gov_kind_vocab,
	      const VocabF& rel_kind_vocab,
	      const VocabT& gov_vocab,
	      const VocabT& rel_vocab);
    template <typename CorpusT>
    void init(int num_templates, int num_slots,
	      CorpusT* corpus,
	      const VocabF& gov_kind_vocab,
	      const VocabF& rel_kind_vocab,
	      const VocabT& gov_vocab,
	      const VocabT& rel_vocab);
    // template <typename G, typename R> void heldout_init(int num_templates, 
    // 							const std::vector<int>& num_slots,
    // 							const Vocabulary<G>& gov_kind_vocab,
    // 							const Vocabulary<R>& rel_kind_vocab);
    int sample_template(const int curr_slot_val);
    int sample_slot(const int curr_template_val);
    int sample_gov_kind(const int curr_template_val);
    int sample_rel_kind(const int curr_template_val, const int curr_slot_val);
    void adjust_for_slot(int templ_val, int slot_val, int di,
			 int ei, int by);
    void adjust_for_template(int val, int di, int ei, int by);
    void unassign_template(int prev, int di, int ei);
    void assign_template(int sampled, int di, int ei);
    void unassign_slot(int templ, int prev, int di, int ei);
    void assign_slot(int templ, int sampled, int di, int ei);
    void adjust_for_gov_kind(int val, int di, int ei, int mi, int by);
    void unassign_gov_kind(int prev, int di, int ei, int mi);
    void assign_gov_kind(int sampled, int di, int ei, int mi);
    void adjust_for_rel_kind(int val, int di, int ei, int mi, int by);
    void unassign_rel_kind(int prev, int di, int ei, int mi);
    void assign_rel_kind(int sampled, int di, int ei, int mi);

    void update_floating_hists(const minsky::Entity& entity,
			       int doc_index, int ent_index);

    template <typename CorpusT, template <typename> class TSW, typename... TSWArgs >
    void learn
    (
     CorpusT* corpus,
     int epoch,
     const StringDiscreteKindPrinter& print_struct,
     DKVWriters* sw_wrapper,
     bool heldout,
     TSWArgs... tsw_args
     );

    const std::vector<std::vector<double> >& doc_template_params() const;
    void sampling_strategy(WithKindsSamplingStrategy* ss);
    template <typename CorpusT>
    void print_labeling(CorpusT* corpus, const std::string& f_name,
			const VocabF& gkv, const VocabF& rkv,
			const VocabT& gv, const VocabT& rv);
    minsky::residual::ResidualUniqueSlots create_minsky_view(const VocabF& gkv_vocab, const VocabF& rkv_vocab, const VocabT& gov_vocab, const VocabT& rel_vocab) const;
    template <typename Protocol>
    void write_usage_posterior(ferrum::thrift::ThriftSmartWriter<Protocol>* tsw,
			       const std::vector<double>& usage);
  private:
    int num_docs_;
    WithKindsSamplingStrategy* sample_strategy_;
    SymmetricHyperparams shp_;

    int num_templates_;
    std::vector<int> num_slots_;
    int num_gov_kinds_;
    int num_rel_kinds_;

    // assignment arrays
    std::vector< std::vector<int> > template_assignments_;
    std::vector< std::vector<int> > slot_assignments_;
    // NOTE: these may only be trivialized initialized
    // initialization only happens if latent_kinds_ == true
    std::vector< std::vector< std::vector<int> > > gov_kind_assignments_;
    std::vector< std::vector< std::vector<int> > > rel_kind_assignments_;

    typedef std::vector<int>* floating_counter_ptr;
    // count tables
    // num_docs * num_templates: how many times a given template is used in each document
    std::vector<std::vector<int> > c_template_;
    floating_counter_ptr c_doc_template_ptr_;
    // num_templates * num_gov_kinds: how many times a given template is used for a governor
    std::vector<std::vector<int> > c_template_gov_kind_;
    // num_templates: how many times a template occurs
    std::vector<int>  c_template_sum_;
    // num_templates * num_slots_per_template: how many times a slot is used (per template) and by doc
    std::vector< std::vector<int> > c_slot_;
    // num_templates * num_slots_per_template * num rel_kinds
    std::vector<std::vector< std::vector<int> > > c_slot_rel_kind_;

    // and now for the "true" observations
    // 
    std::vector< std::vector<int> > c_gov_kind_gov_;
    std::vector< int > c_gov_kind_gov_sum_;
    std::vector< std::vector<int> > c_rel_kind_rel_;
    std::vector< int > c_rel_kind_rel_sum_;

    gsl_rng_type *which_gsl_rng;
    gsl_rng *rnd_gen_;
    bool own_rng_;

    std::map<int, int> floating_gov_hist_;
    std::map<int, int> floating_rel_hist_;

    dmc::gdmc t_usage_dmc_;
    dmc::mfgdmc rel_kind_dmc_;
    dmc::gdmc gov_kind_dmc_;
    dmc::gdmc slot_dmc_;

    dmc::gdmc gov_kind_gov_dmc_;
    dmc::gdmc rel_kind_rel_dmc_;
    
    std::pair<int, int> floating_mention_obs_;

    bool latent_kinds_;

    bool use_lexical_;
    bool use_lexical_set_;

    struct reestimated {
      bool gk, rk, g, r, tu, su;
    };
    void attempt_reestimation(int iteration, reestimated& reest);
    void print_dist_diag(const ferrum::StringDiscreteKindPrinter& print_struct, size_t K);
  };
}

#include "upf/upf_sampling.tcc"

#endif // FERRUM_LIBNAR_CRTLDA_SAMPLING_H_