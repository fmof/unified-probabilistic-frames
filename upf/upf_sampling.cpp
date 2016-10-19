#include "ferrum/crtlda_sampling.hpp"
#include "ferrum/version.hpp"

#include "upf/upf_sampling.hpp"
#include "upf/upf_sampling.tcc"

#include <Eigen/Dense>

namespace ferrum {
  CollapsedGibbsDMC::CollapsedGibbsDMC(int nt, int ns, const SymmetricHyperparams& hp) :
    shp_(hp),
    num_templates_(nt), num_slots_(std::vector<int>(nt, ns)),
    which_gsl_rng((gsl_rng_type*)gsl_rng_mt19937),
    rnd_gen_(gsl_rng_alloc((const gsl_rng_type*)which_gsl_rng)),
    own_rng_(true) {
  }
  CollapsedGibbsDMC::~CollapsedGibbsDMC() {
    if(own_rng_) {
      gsl_rng_free(rnd_gen_);
    }
  }

  void CollapsedGibbsDMC::use_lexical(bool b) {
    use_lexical_ = b;
    use_lexical_set_ = true;
  }
  bool CollapsedGibbsDMC::use_lexical() {
    if(!use_lexical_set_) {
      ERROR << "use_lexical is NOT set: this result may be meaningless";
    }
    return use_lexical_;
  }
  
  int CollapsedGibbsDMC::sample_template(const int curr_slot_val) {
    std::vector<double> lps(num_templates_);
    for(int template_idx = 0; template_idx < num_templates_; template_idx++) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = t_usage_dmc_.log_u_conditional(template_idx, *c_doc_template_ptr_,
						    num_docs_, 1);
      for(auto ghiter : floating_gov_hist_) {
	const int gov_index = ghiter.first;
	const int gov_count = ghiter.second;
	inner += gov_dmc_.log_u_conditional(gov_index, c_template_gov_[template_idx],
					    c_template_sum_[template_idx], gov_count);
      }
      inner += slot_dmc_.log_u_conditional(curr_slot_val, c_slot_[template_idx],
					   c_template_sum_[template_idx], 1);
      lps[template_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  int CollapsedGibbsDMC::sample_slot(const int curr_template_val) {
    int num_slots = num_slots_[curr_template_val];
    std::vector<double> lps(num_slots);
    const std::vector<int>& memoized_template_counts = c_slot_[curr_template_val];
    const int memoized_template_sum = c_template_sum_[curr_template_val];
    for(int slot_idx = 0; slot_idx < num_slots; ++slot_idx) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = slot_dmc_.log_u_conditional(slot_idx, memoized_template_counts,
						 memoized_template_sum, 1);
      for(auto rhiter : floating_rel_hist_) {
	const int rel_index = rhiter.first;
	const int rel_count = rhiter.second;
	inner += rel_dmc_.log_u_conditional(rel_index, c_slot_rel_[curr_template_val][slot_idx],
					    c_slot_[curr_template_val][slot_idx], rel_count);
      }
      lps[slot_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  void CollapsedGibbsDMC::adjust_for_slot(int templ_val, int slot_val, int di,
					  int ei, int by) {
    c_slot_[templ_val][slot_val] += by;
    for(auto rhiter : floating_rel_hist_) {
      c_slot_rel_[templ_val][slot_val][rhiter.first] += by*rhiter.second;
    }
    //c_slot_sum_[templ_val][slot_val] += by;
  }
  void CollapsedGibbsDMC::adjust_for_template(int val, int di, int ei, int by) {
    const int slot_assign = slot_assignments_[di][ei];
    c_template_[di][val] += by;
    c_template_sum_[val] += by;
    // iterate through floating gov hist
    for(auto ghiter : floating_gov_hist_) {
      c_template_gov_[val][ghiter.first] += by*ghiter.second;
    }
    adjust_for_slot(val, slot_assign, di, ei, by);
  }
  void CollapsedGibbsDMC::unassign_template(int prev, int di, int ei) {
    adjust_for_template(prev, di, ei, -1);
  }
  void CollapsedGibbsDMC::assign_template(int sampled, int di, int ei) {
    adjust_for_template(sampled, di, ei, 1);
  }
  void CollapsedGibbsDMC::unassign_slot(int templ, int prev, int di, int ei) {
    adjust_for_slot(templ, prev, di, ei, -1);
  }
  void CollapsedGibbsDMC::assign_slot(int templ, int sampled, int di, int ei) {
    adjust_for_slot(templ, sampled, di, ei, 1);
  }

  void CollapsedGibbsDMC::update_floating_hists(const minsky::Entity& entity) {
    floating_gov_hist_ = minsky::predicate_histogram_on_level(entity, annot_level_);
    floating_rel_hist_ = minsky::relation_histogram_on_level(entity, annot_level_);
  }

  std::vector< std::vector<double> > CollapsedGibbsDMC::doc_template_params() {
    return t_usage_dmc_.collapsed_params();
  }

  //setters
  void CollapsedGibbsDMC::sampling_strategy(SamplingStrategy* ss) {
    sample_strategy_ = ss;
  }

  void
  CollapsedGibbsDMC::print_dist_diag(const ferrum::StringDiscreteKindPrinter& print_struct,
				     size_t K) {
    Eigen::VectorXd entropies(num_docs_);
    int i = 0;
    std::vector<double> avg_use(num_templates_, 0.0);
    for(const auto& usage : t_usage_dmc_.collapsed_params()) {
      double norm = ferrum::sum(usage);
      std::vector<double> prop = ferrum::scalar_product(1.0/norm, usage);
      entropies(i++) = dmc::cat::entropy(prop);
      ferrum::sum_in_first(&avg_use, prop);
    }
    ferrum::scalar_product(1.0/(double)num_docs_, &avg_use);
    INFO << "Template usage summary for " << num_docs_ << " documents: EntropyMean = " << entropies.mean() << ", EntropyVariance = " << ferrum::variance(entropies);
    std::stringstream ss;
    ferrum::to_stringstream(avg_use, ss);
    INFO << "Average p(template) = (" << ss.str() << ")";
    // print the topics now
    std::vector<std::vector<double> > reest_slot_usage = slot_dmc_.collapsed_params();
    int topic_id;
    {
      topic_id = 0;
      for(const auto& topic : reest_slot_usage) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = topic.size();
	std::stringstream stream;
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << which << " (" << topic[which]/norm << ") ";
	}
	INFO << "Template/Slot " << (topic_id++) << ": " << stream.str();
      }
    }
    if(print_struct.g_vocab != NULL) {
      auto* vocab_ptr = print_struct.g_vocab;
      topic_id = 0;
      for(const auto& topic : gov_dmc_.collapsed_params()) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = K > topic.size() ? topic.size() : K;
	std::stringstream stream;
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << vocab_ptr->word(which) << " (" << topic[which]/norm << ") ";
	}
	INFO << "Template/Verb " << (topic_id++) << ": " << stream.str();
      }
    }
    if(print_struct.r_vocab != NULL) {
      auto* vocab_ptr = print_struct.r_vocab;
      int t_id = 0;
      auto rel_p = rel_dmc_.collapsed_params(); // T x S x R
      for(t_id = 0; t_id < num_templates_; ++t_id) {
	const auto& s_topic = reest_slot_usage[t_id];
	std::vector<size_t> sorted_s_topic = ferrum::sort_indices(s_topic, false);
	const auto& s_r_p = rel_p[t_id];
	int s_id = 0;
	for(size_t j = 0; j < sorted_s_topic.size(); ++j) { // iterate number of slot
	  const size_t s_which = sorted_s_topic[j];
	  const auto& topic = s_r_p[ s_which ];
	  std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	  size_t max_K = K > topic.size() ? topic.size() : K;
	  std::stringstream stream;
	  for(size_t i = 0; i < max_K; ++i) {
	    size_t which = sorted_topic[i];
	    stream << vocab_ptr->word(which) << " (" << topic[which] << ") ";
	  }
	  INFO << "Template " << t_id << ", Slot/Rel " << (s_which) << " (" << reest_slot_usage[t_id][s_which]  << "): " << stream.str();
	  ++s_id;
	}
      }
    }
  }


  minsky::residual::ResidualUniqueSlots
  CollapsedGibbsDMC::create_minsky_view
  (
   const VocabT& gov_vocab,
   const VocabT& rel_vocab
   ) const {
    using namespace minsky::residual;
    using namespace minsky;
    ResidualUniqueSlots rtm;
    minsky::Vocab mgv = gov_vocab.minskify();
    rtm.vocabularies.push_back(mgv);
    minsky::Vocab mgr = rel_vocab.minskify();
    rtm.vocabularies.push_back(mgr);
    rtm.__set_predicate_hyper(gov_dmc_.hyperparameters());
    rtm.__set_relation_hyper(rel_dmc_.hyperparameters());    // fix
    rtm.__set_slot_hyper(slot_dmc_.hyperparameters()); // fix
    const auto& gov_params = gov_dmc_.collapsed_params(); // T x V
    const auto& slot_use = slot_dmc_.collapsed_params(); // T x S
    const auto& rel_params = rel_dmc_.collapsed_params(); // T x S x R
    for(size_t ti = 0; ti < (size_t)num_templates_; ++ti) {
      SituationTemplate st;
      {
	minsky::Frame f;
	minsky::Distribution d;
	d.__set_support_size(gov_vocab.num_words());
	minsky::Weights w;
	w.__set_normalized(gov_params[ti]);
	d.__set_weights(w);
	d.__set_vocab_idx(0);
	f.__set_distr(d);
	st.__set_predicate_frame(f);
      }
      // now handle the slots
      {
	SituationSlotFrame ssf;
	// per-template usage
	{
	  Distribution d;
	  d.__set_support_size(num_slots_[ti]);
	  Weights w;
	  w.__set_normalized(slot_use[ti]);
	  d.__set_weights(w);
	  ssf.__set_usage(d);
	}
	// role probs
	{
	  for(size_t si = 0; si < (size_t)num_slots_[ti]; ++si) {
	    SituationSlot ss;
	    Frame sr_frame;
	    Distribution d;
	    d.__set_support_size(rel_vocab.num_words());
	    Weights w;
	    w.__set_normalized(rel_params[ti][si]);
	    d.__set_weights(w);
	    d.__set_vocab_idx(1);
	    sr_frame.__set_distr(d);
	    ss.__set_role(sr_frame);
	    ssf.slots.push_back(ss);
	    ssf.__isset.slots = true;
	  }
	}
	st.__set_unique_slots(ssf);
      }      
      rtm.thematics.push_back(st);
    }
    return rtm;
  }

  template <typename P>
  void CollapsedGibbsDMC::write_usage_posterior
  (
   ferrum::thrift::ThriftSmartWriter<P>* tsw,
   const std::vector<double>& usage_post
   ) {
    concrete::CommunicationTagging ct;
    concrete::util::uuid_factory uf;
    // required: uuid
    concrete::util::add_uuid(ct, uf);
    // required: metadata
    concrete::util::add_metadata(ct,
				 "UPF-syntax-only template induction: " + ferrum::FERRUM_GIT_SHA);
    ct.__set_taggingType("upf-syntax-templates");
    // loop through a (sorted) posterior, setting both the tagList and confidenceList
    const double norm = ferrum::sum(usage_post);
    std::vector<size_t> sorted_topic = ferrum::sort_indices(usage_post, false);
    ct.__isset.tagList = true;
    ct.__isset.confidenceList = true;
    ct.tagList.resize(num_templates_);
    ct.confidenceList.resize(num_templates_);
    std::vector<std::string>& tl = ct.tagList;
    std::vector<double>& cl = ct.confidenceList;
    for(size_t i = 0; i < (size_t)num_templates_; ++i) {
      size_t which = sorted_topic[i];
      tl[i] = std::to_string(which);
      cl[i] = (double)usage_post[which] / norm ;
    }
    ferrum::thrift::save<P, concrete::CommunicationTagging>(tsw, ct);
  }

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  CollapsedGibbsDMCWithKinds::CollapsedGibbsDMCWithKinds
  (int nt, int ns, int nf, int nr, const SymmetricHyperparams& hp,
   bool latent_kinds) :
    shp_(hp),
    num_templates_(nt),
    num_slots_(std::vector<int>(nt, ns)),
    num_gov_kinds_(nf),
    num_rel_kinds_(nr),
    which_gsl_rng((gsl_rng_type*)gsl_rng_mt19937),
    rnd_gen_(gsl_rng_alloc((const gsl_rng_type*)which_gsl_rng)),
    own_rng_(true),
    latent_kinds_(latent_kinds) {
  }
  CollapsedGibbsDMCWithKinds::~CollapsedGibbsDMCWithKinds() {
    if(own_rng_) {
      gsl_rng_free(rnd_gen_);
    }
  }

  void CollapsedGibbsDMCWithKinds::rnd_gen(gsl_rng* rng) {
    if(own_rng_) {
      gsl_rng_free(rnd_gen_);
    }
    own_rng_ = false;
    rnd_gen_ = rng;
  }

  void CollapsedGibbsDMCWithKinds::latent_kinds(bool b) {
    latent_kinds_ = b;
  }
  void CollapsedGibbsDMCWithKinds::use_lexical(bool b) {
    use_lexical_ = b;
    use_lexical_set_ = true;
  }
  bool CollapsedGibbsDMCWithKinds::use_lexical() {
    if(!use_lexical_set_) {
      ERROR << "use_lexical is NOT set: this result may be meaningless";
    }
    return use_lexical_;
  }

  int CollapsedGibbsDMCWithKinds::sample_template(const int curr_slot_val) {
    std::vector<double> lps(num_templates_);
    for(int template_idx = 0; template_idx < num_templates_; template_idx++) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = t_usage_dmc_.log_u_conditional(template_idx, *c_doc_template_ptr_,
						    num_docs_, 1);
      for(auto ghiter : floating_gov_hist_) {
	const int gov_index = ghiter.first;
	const int gov_count = ghiter.second;
	inner += gov_kind_dmc_.log_u_conditional(gov_index, c_template_gov_kind_[template_idx],
						 c_template_sum_[template_idx], gov_count);
      }
      inner += slot_dmc_.log_u_conditional(curr_slot_val, c_slot_[template_idx],
					   c_template_sum_[template_idx], 1);
      lps[template_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  int CollapsedGibbsDMCWithKinds::sample_slot(const int curr_template_val) {
    int num_slots = num_slots_[curr_template_val];
    std::vector<double> lps(num_slots);
    const std::vector<int>& memoized_template_counts = c_slot_[curr_template_val];
    const int memoized_template_sum = c_template_sum_[curr_template_val];
    for(int slot_idx = 0; slot_idx < num_slots; ++slot_idx) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = slot_dmc_.log_u_conditional(slot_idx, memoized_template_counts,
						 memoized_template_sum, 1);
      for(auto rhiter : floating_rel_hist_) {
	const int rel_index = rhiter.first;
	const int rel_count = rhiter.second;
	inner += rel_kind_dmc_.log_u_conditional(rel_index, c_slot_rel_kind_[curr_template_val][slot_idx],
						 c_slot_[curr_template_val][slot_idx], rel_count);
      }
      lps[slot_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  int CollapsedGibbsDMCWithKinds::sample_gov_kind(const int curr_template_val) {
    std::vector<double> lps(num_gov_kinds_);
    floating_counter_ptr c_t_gk = &(c_template_gov_kind_[curr_template_val]);
    const int template_marginal = c_template_sum_[curr_template_val];
    for(int gk_idx = 0; gk_idx < num_gov_kinds_; gk_idx++) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = gov_kind_dmc_.log_u_conditional(gk_idx, *c_t_gk, template_marginal, 1);
      inner += gov_kind_gov_dmc_.log_u_conditional(floating_mention_obs_.first, c_gov_kind_gov_[gk_idx],
						   c_gov_kind_gov_sum_[gk_idx], 1);
      lps[gk_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  int CollapsedGibbsDMCWithKinds::sample_rel_kind(const int curr_template_val, const int curr_slot_val) {
    std::vector<double> lps(num_rel_kinds_);
    floating_counter_ptr c_t_rk = &(c_slot_rel_kind_[curr_template_val][curr_slot_val]);
    const int slot_marginal = c_slot_[curr_template_val][curr_slot_val];
    for(int rk_idx = 0; rk_idx < num_rel_kinds_; rk_idx++) {
      // inner is a log-prob, so addition **here** is multiplication **there**
      double inner = rel_kind_dmc_.log_u_conditional(rk_idx, *c_t_rk, slot_marginal, 1);
      inner += rel_kind_rel_dmc_.log_u_conditional(floating_mention_obs_.second, c_rel_kind_rel_[rk_idx],
						   c_rel_kind_rel_sum_[rk_idx], 1);
      lps[rk_idx] = inner;
    }
    return dmc::cat::log_u_sample(lps);
  }
  void CollapsedGibbsDMCWithKinds::adjust_for_slot(int templ_val, int slot_val, int di,
						   int ei, int by) {
    c_slot_[templ_val][slot_val] += by;
    for(auto rhiter : floating_rel_hist_) {
      c_slot_rel_kind_[templ_val][slot_val][rhiter.first] += by*rhiter.second;
    }
  }
  void CollapsedGibbsDMCWithKinds::adjust_for_template(int val, int di, int ei, int by) {
    const int slot_assign = slot_assignments_[di][ei];
    c_template_[di][val] += by;
    c_template_sum_[val] += by;
    // iterate through floating gov hist
    for(auto ghiter : floating_gov_hist_) {
      c_template_gov_kind_[val][ghiter.first] += by*ghiter.second;
    }
    adjust_for_slot(val, slot_assign, di, ei, by);
  }
  void CollapsedGibbsDMCWithKinds::unassign_template(int prev, int di, int ei) {
    adjust_for_template(prev, di, ei, -1);
  }
  void CollapsedGibbsDMCWithKinds::assign_template(int sampled, int di, int ei) {
    adjust_for_template(sampled, di, ei, 1);
  }
  void CollapsedGibbsDMCWithKinds::unassign_slot(int templ, int prev, int di, int ei) {
    adjust_for_slot(templ, prev, di, ei, -1);
  }
  void CollapsedGibbsDMCWithKinds::assign_slot(int templ, int sampled, int di, int ei) {
    adjust_for_slot(templ, sampled, di, ei, 1);
  }
  void CollapsedGibbsDMCWithKinds::adjust_for_gov_kind(int val, int di, int ei, int mi, int by) {
    const int templ_assign = template_assignments_[di][ei];
    c_template_gov_kind_[templ_assign][val] += by;
    c_gov_kind_gov_sum_[val] += by;
    c_gov_kind_gov_[val][floating_mention_obs_.first] += by;
    // NOTE: It would be nice to say that the relation kind depends on the
    // governor kind. However, that introduces **a lot** of parameters into
    // the model, and I worry about over-fitting. This does mean that for a 
    // gov kind of FOO, a relation kind of BAR-BAZ could theoretically be 
    // generated.
    // // adjust_for_rel_kind(val, slot_assign, di, ei, by);
  }
  void CollapsedGibbsDMCWithKinds::unassign_gov_kind(int prev, int di, int ei, int mi) {
    adjust_for_gov_kind(prev, di, ei, mi, -1);
  }
  void CollapsedGibbsDMCWithKinds::assign_gov_kind(int sampled, int di, int ei, int mi) {
    adjust_for_gov_kind(sampled, di, ei, mi, 1);
  }
  void CollapsedGibbsDMCWithKinds::adjust_for_rel_kind(int val, int di, int ei, int mi, int by) {
    const int templ_assign = template_assignments_[di][ei];
    const int slot_assign = slot_assignments_[di][ei];
    c_slot_rel_kind_[templ_assign][slot_assign][val] += by;
    c_rel_kind_rel_sum_[val] += by;
    c_rel_kind_rel_[val][floating_mention_obs_.second] += by;
  }
  void CollapsedGibbsDMCWithKinds::unassign_rel_kind(int prev, int di, int ei, int mi) {
    adjust_for_rel_kind(prev, di, ei, mi, -1);
  }
  void CollapsedGibbsDMCWithKinds::assign_rel_kind(int sampled, int di, int ei, int mi) {
    adjust_for_rel_kind(sampled, di, ei, mi, 1);
  }
  void CollapsedGibbsDMCWithKinds::update_floating_hists(const minsky::Entity& entity,
							 int doc_index, int ent_index) {
    // if we're not dealing with latent kinds, then we can just read
    // this off of the entity
    if(!latent_kinds_) {
      floating_gov_hist_ = minsky::predicate_histogram_on_level(entity, minsky::AnnotationLevel::SEMANTIC);
      floating_rel_hist_ = minsky::relation_histogram_on_level(entity, minsky::AnnotationLevel::SEMANTIC);
    } else { //otherwise, we need to construct the histogram ourselves
      const int num_m = entity.mentions.size();
      floating_gov_hist_.clear();
      floating_rel_hist_.clear();
      for(int mi = 0; mi < num_m; ++mi) {
	++floating_gov_hist_[gov_kind_assignments_[doc_index][ent_index][mi]];
	++floating_rel_hist_[rel_kind_assignments_[doc_index][ent_index][mi]];
      }
    }
  }
  void CollapsedGibbsDMCWithKinds::attempt_reestimation(int iteration, reestimated& reest) {
    if(sample_strategy_->reestimate_gov_kind(iteration)) {
      gov_kind_dmc_.reestimate_collapsed_parameters(c_template_gov_kind_);
      reest.gk = true;
    }
    if(sample_strategy_->reestimate_rel_kind(iteration)) {
      rel_kind_dmc_.reestimate_collapsed_parameters(c_slot_rel_kind_);
      reest.rk = true;
    }
    if(sample_strategy_->reestimate_template_usage(iteration)) {
      t_usage_dmc_.reestimate_collapsed_parameters(c_template_);
      reest.tu = true;
    }
    if(sample_strategy_->reestimate_slot_usage(iteration)) {
      slot_dmc_.reestimate_collapsed_parameters(c_slot_);
      reest.su = true;
    }
    if(sample_strategy_->reestimate_gov(iteration)) {
      gov_kind_gov_dmc_.reestimate_collapsed_parameters(c_gov_kind_gov_);
      reest.g = true;
    }
    if(sample_strategy_->reestimate_rel(iteration)) {
      rel_kind_rel_dmc_.reestimate_collapsed_parameters(c_rel_kind_rel_);
      reest.r = true;
    }
  }

  const std::vector<std::vector<double> >&
  CollapsedGibbsDMCWithKinds::doc_template_params() const {
    return t_usage_dmc_.collapsed_params();
  }

  //setters
  void CollapsedGibbsDMCWithKinds::sampling_strategy(WithKindsSamplingStrategy* ss) {
    sample_strategy_ = ss;
  }

  void CollapsedGibbsDMCWithKinds::print_dist_diag
  (const ferrum::StringDiscreteKindPrinter& print_struct,
   size_t K) {
    Eigen::VectorXd entropies(num_docs_);
    int i = 0;
    std::vector<double> avg_use(num_templates_, 0.0);
    for(const auto& usage : t_usage_dmc_.collapsed_params()) {
      double norm = ferrum::sum(usage);
      std::vector<double> prop = ferrum::scalar_product(1.0/norm, usage);
      entropies(i++) = dmc::cat::entropy(prop);
      ferrum::sum_in_first(&avg_use, prop);
    }
    ferrum::scalar_product(1.0/(double)num_docs_, &avg_use);
    INFO << "Template usage summary for " << num_docs_ << " documents: EntropyMean = " << entropies.mean() << ", EntropyVariance = " << ferrum::variance(entropies);
    std::stringstream ss;
    ferrum::to_stringstream(avg_use, ss);
    INFO << "Average p(template) = (" << ss.str() << ")";
    // print the topics now
    std::vector<std::vector<double> > reest_slot_usage = slot_dmc_.collapsed_params();
    int topic_id;
    {
      topic_id = 0;
      for(const auto& topic : reest_slot_usage) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = topic.size();
	std::stringstream stream;
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << which << " (" << topic[which]/norm << ") ";
	}
	INFO << "Template/Slot " << (topic_id++) << ": " << stream.str();
      }
    }
    if(print_struct.gk_vocab != NULL) {
      auto* vocab_ptr = print_struct.gk_vocab;
      topic_id = 0;
      for(const auto& topic : gov_kind_dmc_.collapsed_params()) {
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = K > topic.size() ? topic.size() : K;
	std::stringstream stream;
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << vocab_ptr->word(which) << " (" << topic[which]/norm << ") ";
	}
	INFO << "Template/Semantic Gov " << (topic_id++) << ": " << stream.str();
      }
    }
    if(print_struct.g_vocab != NULL) {
      auto* vocab_ptr = print_struct.g_vocab;
      for(topic_id = 0; topic_id < num_templates_; ++topic_id) {
	const auto& ft = gov_kind_dmc_.collapsed_params(topic_id);
	const size_t s = vocab_ptr->num_words();
	std::vector<double> topic(s, 0.0);
	const auto& gpar = gov_kind_gov_dmc_.collapsed_params();
	for(size_t vi = 0; vi < s; ++vi) {
	  for(size_t fi = 0; fi < ft.size(); ++fi) {
	    topic[vi] += ft[fi] * gpar[fi][vi];
	  }
	}
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = K > topic.size() ? topic.size() : K;
	std::stringstream stream;
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << vocab_ptr->word(which) << " (" << topic[which]/norm << ") ";
	}
	INFO << "Marginalized Template/Verb " << topic_id << ": " << stream.str();
      }
      std::vector<double> gk_avg_use(num_gov_kinds_, 0.0);
      for(const auto& topic : gov_kind_dmc_.collapsed_params()) {
	double norm = ferrum::sum(topic);
	std::vector<double> prop = ferrum::scalar_product(1.0/norm, topic);
	ferrum::sum_in_first(&gk_avg_use, prop);
      }
      std::vector<size_t> sorted_gk_avg_use =
	ferrum::sort_indices(gk_avg_use, false);
      size_t max_gk_K = K > (size_t)num_gov_kinds_ ? (size_t)num_gov_kinds_ : K;
      const auto& gkg_params = gov_kind_gov_dmc_.collapsed_params();
      for(topic_id = 0; topic_id < (int)max_gk_K; ++topic_id) {
	const auto& topic = gkg_params[sorted_gk_avg_use[topic_id]];
	std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	double norm = ferrum::sum(topic);
	size_t max_K = K > topic.size() ? topic.size() : K;
	std::stringstream stream;
	stream << "Semantic/Syntactic Gov " << sorted_gk_avg_use[topic_id];
	if(print_struct.gk_vocab != NULL) {
	  stream << "(" << print_struct.gk_vocab->word(sorted_gk_avg_use[topic_id]) << ")";
	}
	stream << " [avg. use: " << gk_avg_use[sorted_gk_avg_use[topic_id]] << "]";
	stream << ": ";
	for(size_t i = 0; i < max_K; ++i) {
	  size_t which = sorted_topic[i];
	  stream << vocab_ptr->word(which) << " (" << topic[which]/norm << ") ";
	}
	INFO << stream.str();
      }
    }
    if(print_struct.rk_vocab != NULL) {
      auto* vocab_ptr = print_struct.rk_vocab;
      int t_id = 0;
      auto rel_p = rel_kind_dmc_.collapsed_params(); // T x S x R
      for(t_id = 0; t_id < num_templates_; ++t_id) {
	const auto& s_topic = reest_slot_usage[t_id];
	std::vector<size_t> sorted_s_topic = ferrum::sort_indices(s_topic, false);
	const auto& s_r_p = rel_p[t_id];
	int s_id = 0;
	for(size_t j = 0; j < sorted_s_topic.size(); ++j) { // iterate number of slot
	  const size_t s_which = sorted_s_topic[j];
	  const auto& topic = s_r_p[ s_which ];
	  std::vector<size_t> sorted_topic = ferrum::sort_indices(topic, false);
	  size_t max_K = K > topic.size() ? topic.size() : K;
	  std::stringstream stream;
	  for(size_t i = 0; i < max_K; ++i) {
	    size_t which = sorted_topic[i];
	    stream << vocab_ptr->word(which) << " (" << topic[which] << ") ";
	  }
	  INFO << "Template " << t_id << ", Slot/Sem Role " << (s_which) << " (" << reest_slot_usage[t_id][s_which]  << "): " << stream.str();
	  ++s_id;
	}
      }
    }
  }
  
  minsky::residual::ResidualUniqueSlots
  CollapsedGibbsDMCWithKinds::create_minsky_view
  (
   const VocabF& gkv,
   const VocabF& rkv,
   const VocabT& gov_vocab,
   const VocabT& rel_vocab
   ) const {
    using namespace minsky::residual;
    using namespace minsky;
    ResidualUniqueSlots rtm;
    rtm.vocabularies.push_back(gkv.minskify());
    rtm.vocabularies.push_back(rkv.minskify());
    rtm.vocabularies.push_back(gov_vocab.minskify());
    rtm.vocabularies.push_back(rel_vocab.minskify());
    rtm.__set_predicate_hyper(gov_kind_gov_dmc_.hyperparameters());
    rtm.__set_relation_hyper(rel_kind_rel_dmc_.hyperparameters());
    rtm.__set_slot_hyper(slot_dmc_.hyperparameters());
    rtm.__set_sem_predicate_hyper(gov_kind_dmc_.hyperparameters());
    rtm.__set_sem_relation_hyper(rel_kind_dmc_.hyperparameters());
    
    const auto& gk_params = gov_kind_dmc_.collapsed_params(); // T x F
    const auto& gov_params = gov_kind_gov_dmc_.collapsed_params(); // F x V
    const auto& slot_use = slot_dmc_.collapsed_params(); // T x S
    const auto& rk_params = rel_kind_dmc_.collapsed_params(); // T x S x R
    const auto& rel_params = rel_kind_rel_dmc_.collapsed_params(); // R x D
    for(size_t ti = 0; ti < (size_t)num_templates_; ++ti) {
      {
	SituationTemplate st;
	{
	  minsky::Frame f;
	  minsky::Distribution d;
	  d.__set_support_size(gkv.num_words());
	  minsky::Weights w;
	  w.__set_normalized(gk_params[ti]);
	  d.__set_weights(w);
	  d.__set_vocab_idx(0);
	  f.__set_distr(d);
	  st.__set_predicate_frame(f);
	}
	// now handle the slots
	{
	  SituationSlotFrame ssf;
	  // per-template usage
	  {
	    Distribution d;
	    d.__set_support_size(num_slots_[ti]);
	    Weights w;
	    w.__set_normalized(slot_use[ti]);
	    d.__set_weights(w);
	    ssf.__set_usage(d);
	  }
	  // role probs
	  {
	    for(size_t si = 0; si < (size_t)num_slots_[ti]; ++si) {
	      SituationSlot ss;
	      Frame sr_frame;
	      Distribution d;
	      d.__set_support_size(rkv.num_words());
	      Weights w;
	      w.__set_normalized(rk_params[ti][si]);
	      d.__set_weights(w);
	      d.__set_vocab_idx(1);
	      sr_frame.__set_distr(d);
	      ss.__set_role(sr_frame);
	      ssf.slots.push_back(ss);
	      ssf.__isset.slots = true;
	    }
	  }
	  st.__set_unique_slots(ssf);
	}      
	rtm.thematics.push_back(st);
      }
    }
    for(size_t fi = 0; fi < (size_t)gkv.num_words(); ++fi) {
      { // fix
	SituationTemplate sem;
	minsky::Frame f;
	minsky::Distribution d;
	d.__set_support_size(gov_vocab.num_words());
	minsky::Weights w;
	w.__set_normalized(gov_params[fi]);
	d.__set_weights(w);
	d.__set_vocab_idx(2);
	f.__set_distr(d);
	sem.__set_predicate_frame(f);
	// now handle the semantic roles
	rtm.semantics.preds.push_back(sem);
	rtm.__isset.semantics = true;
	rtm.semantics.__isset.preds = true;
      }
    }
    for(size_t ri = 0; ri < (size_t)rkv.num_words(); ++ri) {
      SituationSlot ss;
      Frame sr_frame;
      Distribution d;
      d.__set_support_size(rel_vocab.num_words());
      Weights w;
      w.__set_normalized(rel_params[ri]);
      d.__set_weights(w);
      d.__set_vocab_idx(3);
      sr_frame.__set_distr(d);
      ss.__set_role(sr_frame);
      rtm.semantics.roles.push_back(ss);
      rtm.__isset.semantics = true;
      rtm.semantics.__isset.roles = true;
    }
    return rtm;
  }

  template <typename P>
  void CollapsedGibbsDMCWithKinds::write_usage_posterior
  (
   ferrum::thrift::ThriftSmartWriter<P>* tsw,
   const std::vector<double>& usage_post
   ) {
    concrete::CommunicationTagging ct;
    concrete::util::uuid_factory uf;
    // required: uuid
    concrete::util::add_uuid(ct, uf);
    // required: metadata
    concrete::util::add_metadata(ct,
				 "UPF template induction: " + ferrum::FERRUM_GIT_SHA);
    ct.__set_taggingType("upf-syntax-templates");
    // loop through a (sorted) posterior, setting both the tagList and confidenceList
    const double norm = ferrum::sum(usage_post);
    std::vector<size_t> sorted_topic = ferrum::sort_indices(usage_post, false);
    ct.__isset.tagList = true;
    ct.__isset.confidenceList = true;
    ct.tagList.resize(num_templates_);
    ct.confidenceList.resize(num_templates_);
    std::vector<std::string>& tl = ct.tagList;
    std::vector<double>& cl = ct.confidenceList;
    for(size_t i = 0; i < (size_t)num_templates_; ++i) {
      size_t which = sorted_topic[i];
      tl[i] = std::to_string(which);
      cl[i] = (double)usage_post[which] / norm ;
    }
    ferrum::thrift::save<P, concrete::CommunicationTagging>(tsw, ct);
  }
  
} // end namespace ferrum
		 
template
void ferrum::CollapsedGibbsDMC::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TBinaryProtocol>* tsw,
 const std::vector<double>& usage
 );
template
void ferrum::CollapsedGibbsDMC::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TCompactProtocol>* tsw,
 const std::vector<double>& usage
 );
template
void ferrum::CollapsedGibbsDMC::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TJSONProtocol>* tsw,
 const std::vector<double>& usage
 );

template
void ferrum::CollapsedGibbsDMCWithKinds::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TBinaryProtocol>* tsw,
 const std::vector<double>& usage
 );
template
void ferrum::CollapsedGibbsDMCWithKinds::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TCompactProtocol>* tsw,
 const std::vector<double>& usage
 );
template
void ferrum::CollapsedGibbsDMCWithKinds::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TJSONProtocol>* tsw,
 const std::vector<double>& usage
 );
