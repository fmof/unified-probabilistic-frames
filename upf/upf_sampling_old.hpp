  template <typename M, typename D, 
  	    typename GKV, typename RKV>
   class CollapsedGibbsDMCWithKindsLenient {
  private:
    int num_docs_;
    std::vector<int> num_entities_;
    WithKindsSamplingStrategy* sample_strategy_;

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

    std::vector< std::vector< std::vector< bool > > > mention_kinds_latent_;

    // the model itself
    M* model_;
    InMemoryCorpus<D>* corpus_;

    gsl_rng_type *which_gsl_rng;
    gsl_rng *rnd_gen_;

    std::map<int, int> floating_gov_hist_;
    std::map<int, int> floating_rel_hist_;

    dmc::gdmc t_usage_dmc_;
    dmc::mfgdmc rel_kind_dmc_;
    dmc::gdmc gov_kind_dmc_;
    dmc::gdmc slot_dmc_;

    dmc::gdmc gov_kind_gov_dmc_;
    dmc::gdmc rel_kind_rel_dmc_;
    
    GKV *gov_kind_vocab_ptr_;
    RKV *rel_kind_vocab_ptr_;
    Vocabulary<std::string> actual_gov_vocab_;
    Vocabulary<std::string> actual_rel_vocab_;

    std::pair<int, int> floating_mention_obs_;

    bool latent_kinds_;

    bool use_lexical_;
    bool use_lexical_set_;

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar & num_docs_;
      ar & num_entities_;
      ar & num_templates_;
      ar & num_slots_;
      ar & num_gov_kinds_;
      ar & num_rel_kinds_;

      ar & template_assignments_;
      ar & slot_assignments_;
      ar & gov_kind_assignments_;
      ar & rel_kind_assignments_;

      // count tables
      ar & c_template_;
      ar & c_doc_template_ptr_;
      ar &  c_template_gov_kind_;
      ar & c_slot_;
      ar & c_slot_rel_kind_;

      ar & c_gov_kind_gov_;
      ar & c_rel_kind_rel_;

      ar & t_usage_dmc_;
      ar & rel_kind_dmc_;
      ar & gov_kind_dmc_;
      ar & slot_dmc_;
      ar & gov_kind_gov_dmc_;
      ar & rel_kind_rel_dmc_;
      
      ar & mention_kinds_latent_;

      ar & gov_kind_vocab_ptr_;
      ar & rel_kind_vocab_ptr_;
      ar & actual_gov_vocab_;
      ar & actual_rel_vocab_;

      ar & latent_kinds_;

      ar & use_lexical_;
      ar & use_lexical_set_;
    }

  public:
    void latent_kinds(bool b) {
      latent_kinds_ = b;
    }
    void use_lexical(bool b) {
      use_lexical_ = b;
      use_lexical_set_ = true;
    }
    bool use_lexical() {
      if(!use_lexical_set_) {
	ERROR << "use_lexical is NOT set: this result may be meaningless";
      }
      return use_lexical_;
    }
    CollapsedGibbsDMCWithKindsLenient<M, D, GKV, RKV>() : 
      which_gsl_rng((gsl_rng_type*)gsl_rng_mt19937),
      rnd_gen_(gsl_rng_alloc((const gsl_rng_type*)which_gsl_rng)) {
    }
    CollapsedGibbsDMCWithKindsLenient<M, D, GKV, RKV>(M* model, InMemoryCorpus<D>* corpus,
  					       GKV* gov_kind_vocab, RKV* rel_kind_vocab,
  					       const Vocabulary<std::string>& actual_gov_vocab,
  					       const Vocabulary<std::string>& actual_rel_vocab,
					       bool latent_kinds) : 
    num_docs_(corpus->num_docs()),
      num_entities_(corpus->vec_entity_count()),
				    num_templates_(model->num_templates()), num_slots_(model->num_slots()),
      template_assignments_(std::vector<std::vector<int> >(num_docs_)),
      slot_assignments_(std::vector<std::vector<int> >(num_docs_)), 
      model_(model), corpus_(corpus),
      which_gsl_rng((gsl_rng_type*)gsl_rng_mt19937),
      rnd_gen_(gsl_rng_alloc((const gsl_rng_type*)which_gsl_rng)),
      gov_kind_vocab_ptr_(gov_kind_vocab), rel_kind_vocab_ptr_(rel_kind_vocab),
      actual_gov_vocab_(actual_gov_vocab), actual_rel_vocab_(actual_rel_vocab),
      latent_kinds_(latent_kinds) {
    }
    ~CollapsedGibbsDMCWithKindsLenient() {
    }

    void rnd_gen(gsl_rng* rng) {
      rnd_gen_ = rng;
    }

    void model(M* model) {
      model_ = model;
    }
    void corpus(InMemoryCorpus<D>* corpus) {
      corpus_ = corpus;
      if(num_docs_ != corpus_->num_docs()) {
	WARN << "The number of docs has changed from " << num_docs_ << " to " << corpus_->num_docs() << ". Please call heldout_init().";
      }
      num_docs_ = corpus_->num_docs();
      num_entities_ = corpus_->vec_entity_count();
    }

    template <typename G, typename R> void init(int num_templates, 
  						const std::vector<int>& num_slots,
  						const Vocabulary<G>& gov_kind_vocab,
  						const Vocabulary<R>& rel_kind_vocab) {
      gov_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
      rel_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
      mention_kinds_latent_ = std::vector< std::vector< std::vector<bool> > >(num_docs_);

      num_gov_kinds_ = model_->num_gov_kinds();
      num_rel_kinds_ = model_->num_rel_kinds();
      if(num_gov_kinds_ < 0 || num_rel_kinds_ < 0) {
	ERROR << "There was an issue in getting either the number of latent governor kinds ( " << num_gov_kinds_ << ") or the number of latent relation kinds (" << num_rel_kinds_ << "). Please fix.";
	throw 10;
      }

      t_usage_dmc_ = dmc::gdmc(model_->num_templates(), model_->hyper_theta(), num_docs_);
      gov_kind_dmc_ = dmc::gdmc(gov_kind_vocab.num_words(), model_->hyper_gov_kind(), model_->num_templates());
      rel_kind_dmc_ = dmc::mfgdmc(rel_kind_vocab.num_words(), model_->hyper_rel_kind(), model_->num_templates(), num_slots_);
      gov_kind_gov_dmc_ = dmc::gdmc(actual_gov_vocab_.num_words(), model_->hyper_gov(), gov_kind_vocab.num_words());
      rel_kind_rel_dmc_ = dmc::gdmc(actual_rel_vocab_.num_words(), model_->hyper_rel(), rel_kind_vocab.num_words());
      slot_dmc_ = dmc::gdmc(num_slots_[0], model_->hyper_slot(), model_->num_templates());

      c_template_sum_.resize(num_templates, 0);
      for(int gk = 0; gk < gov_kind_vocab.num_words(); ++gk) {
	c_gov_kind_gov_.push_back(std::vector<int>(actual_gov_vocab_.num_words(), 0));
	c_gov_kind_gov_sum_.push_back(0);
      }
      for(int rk = 0; rk < rel_kind_vocab.num_words(); ++rk) {
	c_rel_kind_rel_.push_back(std::vector<int>(actual_rel_vocab_.num_words(), 0));
	c_rel_kind_rel_sum_.push_back(0);
      }
      for(int t = 0; t < num_templates; t++) {
  	c_template_gov_kind_.push_back(std::vector<int>(gov_kind_vocab.num_words(), 0));
  	c_slot_.push_back(std::vector< int>(num_slots[t]));
  	c_slot_rel_kind_.push_back(std::vector< std::vector< int> >(num_slots[t]));
  	num_slots_[t] = num_slots[t];
  	for(int s = 0; s < num_slots_[t]; s++) {
  	  c_slot_rel_kind_[t][s] = std::vector<int>(rel_kind_vocab.num_words());
  	}
      }
      for(int di = 0; di < num_docs_; di++) {
  	c_template_.push_back(std::vector<int>(num_templates, 0));
  	const D& doc = (*corpus_)[di];
  	const int num_ents = num_entities_[di];
  	template_assignments_[di] = std::vector<int>(num_ents, 0);
  	slot_assignments_[di]     = std::vector<int>(num_ents, 0);
	gov_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
	rel_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
	mention_kinds_latent_.push_back(std::vector< std::vector< bool > >(num_ents));
  	for(int ei = 0; ei < num_ents; ei++){
  	  int sampled_template = gsl_rng_uniform_int(rnd_gen_, num_templates);
  	  template_assignments_[di][ei] = sampled_template;
  	  int sampled_slot = gsl_rng_uniform_int(rnd_gen_, num_slots_[sampled_template]);
  	  slot_assignments_[di][ei] = sampled_slot;
  	  // update count tables
  	  // these are the latent/latent count tables
  	  ++c_template_[di][sampled_template];
  	  ++c_template_sum_[sampled_template];
  	  ++c_slot_[sampled_template][sampled_slot];
  	  // now update the latent/possibly-observed (kind) count tables
  	  auto& entity = doc[ei];
	  const int num_ments = entity.num_mentions();
	  gov_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));
	  rel_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));	  
	  mention_kinds_latent_[di].push_back(std::vector< bool >(num_ments, false));
	  for(int mi = 0; mi < num_ments; ++mi) {
	    const auto& mention = entity.mention(mi);
	    const bool is_lat = mention.latent();
	    mention_kinds_latent_[di][ei][mi] = is_lat;
	    int gov_kind = is_lat ? gsl_rng_uniform_int(rnd_gen_, num_gov_kinds_) : gov_kind_vocab.index(mention.gov().view());
	    ++c_template_gov_kind_[sampled_template][gov_kind];
	    ++c_gov_kind_gov_sum_[gov_kind];
	    gov_kind_assignments_[di][ei][mi] = gov_kind;
	    int gv = actual_gov_vocab_.index(mention.gov().lemma());
	    ++c_gov_kind_gov_[gov_kind][gv];
	    ////////////////////////////
	    int rel_kind = is_lat ? gsl_rng_uniform_int(rnd_gen_, num_rel_kinds_) : rel_kind_vocab.index(mention.rel());
	    ++c_slot_rel_kind_[sampled_template][sampled_slot][rel_kind];
	    ++c_rel_kind_rel_sum_[rel_kind];
	    rel_kind_assignments_[di][ei][mi] = rel_kind;
	    int rv = actual_rel_vocab_.index(mention.rel_str());
	    ++c_rel_kind_rel_[rel_kind][rv];
	  }
	}
      }
    }

    template <typename G, typename R> void init(int num_templates, int num_slots,
  						const Vocabulary<G>& gov_kind_vocab,
  						const Vocabulary<R>& rel_kind_vocab) {
      init(num_templates_, std::vector<int>(num_templates_, num_slots), gov_kind_vocab, rel_kind_vocab);
    }

    template <typename G, typename R> void init(const Vocabulary<G>& gov_kind_vocab,
  						const Vocabulary<R>& rel_kind_vocab) {
      init(num_templates_, num_slots_, gov_kind_vocab, rel_kind_vocab);
    }

    template <typename G, typename R> void heldout_init(int num_templates, 
							const std::vector<int>& num_slots,
							const Vocabulary<G>& gov_kind_vocab,
							const Vocabulary<R>& rel_kind_vocab) {
      gov_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
      rel_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
      mention_kinds_latent_ = std::vector< std::vector< std::vector<bool> > >(num_docs_);

      num_gov_kinds_ = model_->num_gov_kinds();
      num_rel_kinds_ = model_->num_rel_kinds();
      if(num_gov_kinds_ < 0 || num_rel_kinds_ < 0) {
	ERROR << "There was an issue in getting either the number of latent governor kinds ( " << num_gov_kinds_ << ") or the number of latent relation kinds (" << num_rel_kinds_ << "). Please fix.";
	throw 10;
      }

      t_usage_dmc_ = dmc::gdmc(model_->num_templates(), model_->hyper_theta(), num_docs_);
      
      c_gov_kind_gov_.clear();
      c_gov_kind_gov_sum_.clear();
      c_rel_kind_rel_.clear();
      c_rel_kind_rel_sum_.clear();
      c_template_gov_kind_.clear();
      c_slot_.clear();
      c_slot_rel_kind_.clear();
      for(int gk = 0; gk < gov_kind_vocab.num_words(); ++gk) {
      	c_gov_kind_gov_.push_back(std::vector<int>(actual_gov_vocab_.num_words(), 0));
      	c_gov_kind_gov_sum_.push_back(0);
      }
      for(int rk = 0; rk < rel_kind_vocab.num_words(); ++rk) {
      	c_rel_kind_rel_.push_back(std::vector<int>(actual_rel_vocab_.num_words(), 0));
      	c_rel_kind_rel_sum_.push_back(0);
      }
      for(int t = 0; t < num_templates; t++) {
      	c_template_gov_kind_.push_back(std::vector<int>(gov_kind_vocab.num_words(), 0));
      	c_slot_.push_back(std::vector< int>(num_slots[t]));
      	c_slot_rel_kind_.push_back(std::vector< std::vector< int> >(num_slots[t]));
      	num_slots_[t] = num_slots[t];
      	for(int s = 0; s < num_slots_[t]; s++) {
      	  c_slot_rel_kind_[t][s] = std::vector<int>(rel_kind_vocab.num_words());
      	}
      }
      c_template_.clear();
      c_template_sum_.resize(num_templates, 0);
      template_assignments_.clear();
      slot_assignments_.clear();
      for(int di = 0; di < num_docs_; di++) {
  	c_template_.push_back(std::vector<int>(num_templates, 0));
  	const D& doc = (*corpus_)[di];
  	const int num_ents = num_entities_[di];
  	template_assignments_.push_back(std::vector<int>(num_ents, 0));
  	slot_assignments_.push_back(std::vector<int>(num_ents, 0));
	gov_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
	rel_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
	mention_kinds_latent_.push_back(std::vector< std::vector< bool > >(num_ents));
  	for(int ei = 0; ei < num_ents; ei++){
  	  int sampled_template = gsl_rng_uniform_int(rnd_gen_, num_templates);
  	  template_assignments_[di][ei] = sampled_template;
  	  int sampled_slot = gsl_rng_uniform_int(rnd_gen_, num_slots_[sampled_template]);
  	  slot_assignments_[di][ei] = sampled_slot;
  	  // update count tables
  	  // these are the latent/latent count tables
  	  ++c_template_[di][sampled_template];
  	  ++c_template_sum_[sampled_template];
  	  ++c_slot_[sampled_template][sampled_slot];
  	  // now update the latent/possibly-observed (kind) count tables
  	  auto& entity = doc[ei];
	  const int num_ments = entity.num_mentions();
	  gov_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));
	  rel_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));	  
	  mention_kinds_latent_[di].push_back(std::vector< bool >(num_ments, false));
	  for(int mi = 0; mi < num_ments; ++mi) {
	    const auto& mention = entity.mention(mi);
	    const bool is_lat = mention.latent();
	    mention_kinds_latent_[di][ei][mi] = is_lat;
	    int gov_kind = is_lat ? gsl_rng_uniform_int(rnd_gen_, num_gov_kinds_) : gov_kind_vocab.index(mention.gov().view());
	    ++c_template_gov_kind_[sampled_template][gov_kind];
	    ++c_gov_kind_gov_sum_[gov_kind];
	    gov_kind_assignments_[di][ei][mi] = gov_kind;
	    int gv = actual_gov_vocab_.index(mention.gov().lemma());
	    ++c_gov_kind_gov_[gov_kind][gv];
	    ////////////////////////////
	    int rel_kind = is_lat ? gsl_rng_uniform_int(rnd_gen_, num_rel_kinds_) : rel_kind_vocab.index(mention.rel());
	    ++c_slot_rel_kind_[sampled_template][sampled_slot][rel_kind];
	    ++c_rel_kind_rel_sum_[rel_kind];
	    rel_kind_assignments_[di][ei][mi] = rel_kind;
	    int rv = actual_rel_vocab_.index(mention.rel_str());
	    ++c_rel_kind_rel_[rel_kind][rv];
	  }
	}
      }
    }

    void init() {
      if(rnd_gen_ == NULL) {
	ERROR << "You must set the random number generator by calling rnd_gen(gsl_rng*)";
	return;
      }
      init(num_templates_, num_slots_, *gov_kind_vocab_ptr_, *rel_kind_vocab_ptr_);
    }

    void heldout_init() {
      if(rnd_gen_ == NULL) {
	ERROR << "You must set the random number generator by calling rnd_gen(gsl_rng*)";
	return;
      }
      heldout_init(num_templates_, num_slots_, *gov_kind_vocab_ptr_, *rel_kind_vocab_ptr_);
    }

    int sample_template(const int curr_slot_val) {
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
    int sample_slot(const int curr_template_val) {
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
    int sample_gov_kind(const int curr_template_val) {
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
    int sample_rel_kind(const int curr_template_val, const int curr_slot_val) {
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
    inline void adjust_for_slot(int templ_val, int slot_val, int di,
  				int ei, int by) {
      c_slot_[templ_val][slot_val] += by;
      for(auto rhiter : floating_rel_hist_) {
  	c_slot_rel_kind_[templ_val][slot_val][rhiter.first] += by*rhiter.second;
      }
    }
    inline void adjust_for_template(int val, int di, int ei, int by) {
      const int slot_assign = slot_assignments_[di][ei];
      c_template_[di][val] += by;
      c_template_sum_[val] += by;
      // iterate through floating gov hist
      for(auto ghiter : floating_gov_hist_) {
  	c_template_gov_kind_[val][ghiter.first] += by*ghiter.second;
      }
      adjust_for_slot(val, slot_assign, di, ei, by);
    }
    inline void unassign_template(int prev, int di, int ei) {
      adjust_for_template(prev, di, ei, -1);
    }
    inline void assign_template(int sampled, int di, int ei) {
      adjust_for_template(sampled, di, ei, 1);
    }
    inline void unassign_slot(int templ, int prev, int di, int ei) {
      adjust_for_slot(templ, prev, di, ei, -1);
    }
    inline void assign_slot(int templ, int sampled, int di, int ei) {
      adjust_for_slot(templ, sampled, di, ei, 1);
    }

    inline void adjust_for_gov_kind(int val, int di, int ei, int mi, int by) {
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
    inline void unassign_gov_kind(int prev, int di, int ei, int mi) {
      adjust_for_gov_kind(prev, di, ei, mi, -1);
    }
    inline void assign_gov_kind(int sampled, int di, int ei, int mi) {
      adjust_for_gov_kind(sampled, di, ei, mi, 1);
    }
    inline void adjust_for_rel_kind(int val, int di, int ei, int mi, int by) {
      const int templ_assign = template_assignments_[di][ei];
      const int slot_assign = slot_assignments_[di][ei];
      c_slot_rel_kind_[templ_assign][slot_assign][val] += by;
      c_rel_kind_rel_sum_[val] += by;
      c_rel_kind_rel_[val][floating_mention_obs_.second] += by;
    }
    inline void unassign_rel_kind(int prev, int di, int ei, int mi) {
      adjust_for_rel_kind(prev, di, ei, mi, -1);
    }
    inline void assign_rel_kind(int sampled, int di, int ei, int mi) {
      adjust_for_rel_kind(sampled, di, ei, mi, 1);
    }

    template <typename G, typename R>
    inline void update_floating_hists(int doc_index, int ent_index,
  				      const Vocabulary<G>& gov_kind_vocab,
  				      const Vocabulary<R>& rel_kind_vocab) {
      const D& doc = (*corpus_)[doc_index];
      const auto& entity = doc[ent_index];
      // because some mentions may be latentn, we need to construct the histogram ourselves
      const int num_m = entity.num_mentions();
      floating_gov_hist_.clear();
      floating_rel_hist_.clear();
      for(int mi = 0; mi < num_m; ++mi) {
	++floating_gov_hist_[gov_kind_assignments_[doc_index][ent_index][mi]];
	++floating_rel_hist_[rel_kind_assignments_[doc_index][ent_index][mi]];
      }
    }

    template <typename G, typename R> void learn(const Vocabulary<G>& gov_kind_vocab,
  						 const Vocabulary<R>& rel_kind_vocab,
						 CollapsedGibbsDMCWithKindsLenient<M, D, GKV, RKV>* heldout_sampler,
						 int iter_offset = 0) {
      clock_t start = clock();
      clock_t inter_start = start;
      for(int iteration = iter_offset; iteration < iter_offset + sample_strategy_->num_iterations; ++iteration) {
  	for(int di = 0; di < num_docs_; di++) {
  	  const int num_ents = num_entities_[di];
  	  for(int ei = 0; ei < num_ents; ei++) {
  	    bool floating_hists_updated = false;
  	    if(sample_strategy_->sample_template(iteration, di, ei)) {
  	      update_floating_hists(di, ei, gov_kind_vocab, rel_kind_vocab);
  	      floating_hists_updated = true;
  	      // template sampling
  	      {
  		const int prev_template = template_assignments_[di][ei];
  		c_doc_template_ptr_ = &(c_template_[di]);
  		unassign_template(prev_template, di, ei);
  		int sampled = sample_template(slot_assignments_[di][ei]);
  		assign_template(sampled, di, ei);
  		TRACE << "sampled template:" << sampled;
  		//add counts back
  		if(prev_template != sampled) {
  		  template_assignments_[di][ei] = sampled;
  		}
  	      }
  	    }
  	    if(sample_strategy_->sample_slot(iteration, di, ei)) {
  	      if(!floating_hists_updated) {
  		update_floating_hists(di,ei, gov_kind_vocab, rel_kind_vocab);
  	      }
  	      // slot sampling
  	      {
  		const int prev_slot = slot_assignments_[di][ei];
  		const int templ_assignment = template_assignments_[di][ei];
  		unassign_slot(templ_assignment, prev_slot, di, ei);
  		int sampled = sample_slot(templ_assignment);
  		assign_slot(templ_assignment, sampled, di, ei);
  		TRACE << "sampled slot:" << sampled << " of template " << templ_assignment;
  		//add counts back
  		if(prev_slot != sampled) {
  		  slot_assignments_[di][ei] = sampled;
  		}
  	      }
  	    }
	    // now, if the kinds are latent, then we might need to resample them	    
	    const auto& entity = (*corpus_)[di][ei];
	    const int num_ments = entity.num_mentions();
	    for(int mi = 0; mi < num_ments; ++mi) {
	      const auto& mention = entity.mention(mi);
	      const bool is_lat = mention.latent();
	      floating_mention_obs_.first =  actual_gov_vocab_.index(mention.gov().lemma());
	      floating_mention_obs_.second = actual_rel_vocab_.index(mention.rel_str());
	      if(sample_strategy_->sample_gov_kind(iteration, di, ei, mi, is_lat)) {
		const int prev_gov_kind = gov_kind_assignments_[di][ei][mi];
		unassign_gov_kind(prev_gov_kind, di, ei, mi);
		int sampled = sample_gov_kind(template_assignments_[di][ei]);
		assign_gov_kind(sampled, di, ei, mi);
		TRACE << "sampled gov kind:" << sampled;
		//add counts back
		if(prev_gov_kind != sampled) {
		  gov_kind_assignments_[di][ei][mi] = sampled;
		}
	      }
	      if(sample_strategy_->sample_rel_kind(iteration, di, ei, mi, is_lat)) {
		const int prev_rel_kind = rel_kind_assignments_[di][ei][mi];
		unassign_rel_kind(prev_rel_kind, di, ei, mi);
		int sampled = sample_rel_kind(template_assignments_[di][ei], slot_assignments_[di][ei]);
		assign_rel_kind(sampled, di, ei, mi);
		TRACE << "sampled rel kind:" << sampled;
		//add counts back
		if(prev_rel_kind != sampled) {
		  rel_kind_assignments_[di][ei][mi] = sampled;
		}
	      }
	    }	    
  	  }
  	}
  	clock_t end = clock();
  	double duration = (double)(end - inter_start)/(double)(CLOCKS_PER_SEC);
  	INFO << "Iteration " << (iteration+1) << " took " << duration << " [s]";
  	inter_start = end;
	if(sample_strategy_->reestimate_hyperparameters(iteration)) {
	  CollapsedGibbsDMCWithKindsLenient<M, D, GKV, RKV>* which_sampler = 
	    heldout_sampler != NULL ? heldout_sampler : this;
	  if(heldout_sampler != NULL) {
	    heldout_sampler->learn(gov_kind_vocab, rel_kind_vocab);
	  }
	  t_usage_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_);
	  slot_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_);
	  gov_kind_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_gov_kind_);
	  rel_kind_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_rel_kind_);
	  gov_kind_gov_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_gov_kind_gov_);
	  rel_kind_rel_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_rel_kind_rel_);
	}
	attempt_reestimation(iteration);
      }
    }

    template <typename G, typename R> void learn(const Vocabulary<G>& gov_kind_vocab,
  						 const Vocabulary<R>& rel_kind_vocab,
						 int iter_offset = 0) {
      learn(gov_kind_vocab, rel_kind_vocab, NULL, iter_offset);
    }

    void learn() {
      learn(*gov_kind_vocab_ptr_, *rel_kind_vocab_ptr_, NULL, 0);
    }

    void attempt_reestimation(int iteration) {
      if(sample_strategy_->reestimate_gov_kind(iteration)) {
	gov_kind_dmc_.reestimate_collapsed_parameters(c_template_gov_kind_);
      }
      if(sample_strategy_->reestimate_rel_kind(iteration)) {
	rel_kind_dmc_.reestimate_collapsed_parameters(c_slot_rel_kind_);
      }
      if(sample_strategy_->reestimate_template_usage(iteration)) {
	t_usage_dmc_.reestimate_collapsed_parameters(c_template_);
      }
      if(sample_strategy_->reestimate_slot_usage(iteration)) {
	slot_dmc_.reestimate_collapsed_parameters(c_slot_);
      }
      if(sample_strategy_->reestimate_gov(iteration)) {
	gov_kind_gov_dmc_.reestimate_collapsed_parameters(c_gov_kind_gov_);
      }
      if(sample_strategy_->reestimate_rel(iteration)) {
	rel_kind_rel_dmc_.reestimate_collapsed_parameters(c_rel_kind_rel_);
      }
    }

    M reconstruct_model() {
      M model(this->num_templates_, this->num_slots_);
      model.hyper_theta(t_usage_dmc_.hyperparameters());
      model.hyper_slot(slot_dmc_.hyperparameters());
      model.hyper_gov(gov_kind_gov_dmc_.hyperparameters());
      model.hyper_rel(rel_kind_rel_dmc_.hyperparameters());
      model.hyper_gov_kind(gov_kind_dmc_.hyperparameters());
      model.hyper_rel_kind(rel_kind_dmc_.hyperparameters());
      transfer_learned_parameters(&model);
      return model;
    }

    // transfer the learned parameters back to the model
    void transfer_learned_parameters(M* model) {
      model->prior_template_usage(t_usage_dmc_.collapsed_params());
      model->prior_gov_kind(gov_kind_dmc_.collapsed_params());      
      model->prior_slot_usage(slot_dmc_.collapsed_params());
      model->prior_rel_kind(rel_kind_dmc_.collapsed_params());
      model->prior_gov(gov_kind_gov_dmc_.collapsed_params()); 
      model->prior_rel(rel_kind_rel_dmc_.collapsed_params());
    }
    void transfer_learned_parameters() {
      transfer_learned_parameters(model_);
    }

    GKV gov_kind_vocab() {
      return *gov_kind_vocab_ptr_;
    }
    RKV rel_kind_vocab() {
      return *rel_kind_vocab_ptr_;
    }
    Vocabulary<std::string> gov_vocab() {
      return actual_gov_vocab_;
    }
    Vocabulary<std::string> rel_vocab() {
      return actual_rel_vocab_;
    }

    std::vector<std::vector<double> > doc_template_params() {
      return t_usage_dmc_.collapsed_params();
    }

    //setters
    void sampling_strategy(WithKindsSamplingStrategy* ss) {
      sample_strategy_ = ss;
    }

    void print_labeling(const std::string& f_name) {
      std::ofstream myfile;
      myfile.open(f_name);
      myfile << "doc_id\tent\ttemplate\tslot\tmention\tgov_o\trel_o\tgov_v\trel_v\tgov_k\trel_k\tlatent_kinds\n";
      for(int di = 0; di < num_docs_; di++) {
	const auto& doc = (*corpus_)[di];
	const int num_ents = num_entities_[di];
	for(int ei = 0; ei < num_ents; ei++) {
	  const int prev_template = template_assignments_[di][ei];
	  const int prev_slot = slot_assignments_[di][ei];
	  std::string doc_ent_str = doc.id + "\t" + std::to_string(ei) + 
	    "\t" + std::to_string(prev_template) + "\t" + std::to_string(prev_slot) + "\t";
	  const auto& entity = doc[ei];
	  const int num_ments = entity.num_mentions();
	  for(int mi = 0; mi < num_ments; ++mi) {
	    const auto& mention = entity.mention(mi);
	    std::string actual_gov = mention.gov().lemma();
	    std::string actual_rel = mention.rel_str();
	    myfile << doc_ent_str;
	    myfile << mi << "\t" << actual_gov << "\t" << actual_rel << "\t";
	    myfile << mention.gov().view() << "\t" << mention.rel() << "\t";
	    myfile << gov_kind_assignments_[di][ei][mi] << "\t" << rel_kind_assignments_[di][ei][mi] << "\ttrue\n";
	  }
	}
      }
      myfile.close();
      INFO << "Wrote labeling to " << f_name;
    }
  };
