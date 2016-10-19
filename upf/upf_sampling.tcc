#ifndef UPF_SAMPLING_TCC_
#define UPF_SAMPLING_TCC_

#include "upf/upf_sampling.hpp"
#include "ferrum/timer.hpp"

namespace ferrum {
  template <typename CorpusT>
  void CollapsedGibbsDMC::init(int num_templates,
			       const std::vector<int>& num_slots,
			       CorpusT* corpus,
			       const VocabT& gov_vocab,
			       const VocabT& rel_vocab) {
    num_docs_ = corpus->num_docs();
    num_templates_ = num_templates;
    c_template_sum_.resize(num_templates_, 0);
    template_assignments_ = std::vector<std::vector<int> >(num_docs_);
    slot_assignments_ = std::vector<std::vector<int> >(num_docs_);
    num_slots_ = num_slots;

    t_usage_dmc_ = dmc::gdmc(num_templates_,
			     std::vector<double>(num_templates_, shp_.h_theta),
			     num_docs_);
    slot_dmc_ = dmc::gdmc(num_slots_[0],
			  std::vector<double>(num_slots_[0], shp_.h_slot),
			  num_templates_);
    gov_dmc_ = dmc::gdmc(gov_vocab.num_words(),
			 std::vector<double>(gov_vocab.num_words(), shp_.h_gov),
			 num_templates_);
    rel_dmc_ = dmc::mfgdmc(rel_vocab.num_words(),
			   std::vector<double>(rel_vocab.num_words(), shp_.h_rel),
    			   num_templates_, num_slots_);

    // allocate (and 0-init) global count tables
    for(int t = 0; t < num_templates_; t++) {
      c_template_gov_.push_back(std::vector<int>(gov_vocab.num_words(), 0));
      c_slot_.push_back(std::vector< int>(num_slots_[t]));
      c_slot_rel_.push_back(std::vector< std::vector< int> >(num_slots_[t]));
      for(int s = 0; s < num_slots_[t]; s++) {
	c_slot_rel_[t][s] = std::vector<int>(rel_vocab.num_words());	  
      }
    }
    // allocate (and 0-init) local count tables
    for(int di = 0; di < num_docs_; di++){
      c_template_.push_back(std::vector<int>(num_templates_, 0));
      //const DocT& doc = corpus->operator[](di);
      DocT doc;
      corpus->fill(di, &doc);
      const int num_ents = ferrum::num_entities(doc);
      template_assignments_[di] = std::vector<int>(num_ents, 0);
      slot_assignments_[di]     = std::vector<int>(num_ents, 0);
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
	// now update the latent/observed count tables
	auto& entity = doc.entities[ei];
	std::map<int, int> gov_hist = minsky::predicate_histogram_on_level(entity, annot_level_);
	for(auto ghiter : gov_hist) {
	  c_template_gov_[sampled_template][ghiter.first] += ghiter.second;
	}
	std::map<int, int> rel_hist = minsky::relation_histogram_on_level(entity, annot_level_);
	for(auto rhiter : rel_hist) {
	  c_slot_rel_[sampled_template][sampled_slot][rhiter.first] += rhiter.second;
	}
      }
    }
  }

  template <typename CorpusT>
  void CollapsedGibbsDMC::init(int num_templates, int num_slots,
			       CorpusT* corpus,
			       const VocabT& gov_vocab,
			       const VocabT& rel_vocab) {
    init(num_templates_, std::vector<int>(num_templates_, num_slots), corpus, gov_vocab, rel_vocab);
  }

  template <typename CorpusT>
  void CollapsedGibbsDMC::init(CorpusT* corpus,
			       const VocabT& gov_vocab,
			       const VocabT& rel_vocab) {
    init(num_templates_, num_slots_, corpus, gov_vocab, rel_vocab);
  }

  // template <typename G, typename R> void heldout_init(int num_templates, 
  // 						      const std::vector<int>& num_slots,
  // 						      const Vocabulary<G>& gov_vocab,
  // 						      const Vocabulary<R>& rel_vocab) {
  //   c_template_sum_.resize(num_templates, 0);
  //   c_template_gov_.clear();
  //   c_slot_.clear();
  //   c_slot_rel_.clear();
  //   c_template_.clear();
  //   template_assignments_.clear();
  //   slot_assignments_.clear();
  //   t_usage_dmc_ = dmc::gdmc(model_->num_templates(), model_->hyper_theta(), num_docs_);
  //   for(int t = 0; t < num_templates; t++) {
  //     c_template_gov_.push_back(std::vector<int>(gov_vocab.num_words(), 0));
  //     c_slot_.push_back(std::vector< int>(num_slots[t]));
  //     c_slot_rel_.push_back(std::vector< std::vector< int> >(num_slots[t]));
  //     num_slots_[t] = num_slots[t];
  //     for(int s = 0; s < num_slots_[t]; s++) {
  // 	c_slot_rel_[t][s] = std::vector<int>(rel_vocab.num_words());	  
  //     }
  //   }
  //   for(int di = 0; di < num_docs_; di++){
  //     c_template_.push_back(std::vector<int>(num_templates, 0));
  //     const D& doc = (*corpus_)[di];
  //     const int num_ents = num_entities_[di];
  //     template_assignments_.push_back(std::vector<int>(num_ents, 0));
  //     slot_assignments_.push_back(std::vector<int>(num_ents, 0));
  //     for(int ei = 0; ei < num_ents; ei++){
  // 	int sampled_template = gsl_rng_uniform_int(rnd_gen_, num_templates);
  // 	template_assignments_[di][ei] = sampled_template;
  // 	int sampled_slot = gsl_rng_uniform_int(rnd_gen_, num_slots_[sampled_template]);
  // 	slot_assignments_[di][ei] = sampled_slot;
  // 	// update count tables
  // 	// these are the latent/latent count tables
  // 	++c_template_[di][sampled_template];
  // 	++c_template_sum_[sampled_template];
  // 	++c_slot_[sampled_template][sampled_slot];
  // 	// now update the latent/observed count tables
  // 	auto& entity = doc[ei];
  // 	std::map<G, int> gov_hist = entity.gov_histogram();
  // 	for(auto ghiter : gov_hist) {
  // 	  c_template_gov_[sampled_template][gov_vocab.index(ghiter.first)] += ghiter.second;
  // 	}
  // 	std::map<G, int> rel_hist = entity.rel_histogram();
  // 	for(auto rhiter : rel_hist) {
  // 	  int ri = rel_vocab.index(rhiter.first);
  // 	  c_slot_rel_[sampled_template][sampled_slot][ri] += rhiter.second;
  // 	}
  //     }
  //   }
  // }

  template <typename CorpusT, template <typename> class TSW, typename... TSWArgs >
  void CollapsedGibbsDMC::learn
  (
   CorpusT* corpus,
   int epoch,
   const StringDiscreteKindPrinter& print_struct,
   DKVWriters* sw_wrapper,
   bool heldout,
   TSWArgs... tsw_args
   ) {
    ferrum::Timer timer("Collapsed Gibbs");
    struct {
      bool g, r, tu, su;
    } reestimated;
    reestimated = {heldout, heldout, false, heldout};
    ferrum::ConcreteSituationLabeler<ThriftProtocol, TSW> sit_lab(true);
    sit_lab.template make_with_args< std::string, TSWArgs... >("template_props.tcompact", tsw_args...);
    sit_lab.corpus(corpus);
    for(int iteration = 0; iteration < sample_strategy_->num_iterations; ++iteration) {
      for(int di = 0; di < num_docs_; di++) {
	DocT doc;
	corpus->fill(di, &doc);
	const int num_ents = doc.entities.size();
	for(int ei = 0; ei < num_ents; ei++) {
	  bool floating_hists_updated = false;
	  const minsky::Entity& ent = doc.entities[ei];
	  if(sample_strategy_->sample_template(iteration, di, ei)) {
	    update_floating_hists(ent);
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
	      update_floating_hists(ent);
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
	}
      }
      auto duration = timer.lap();
      INFO << "Iteration " << (iteration+1) << " took " << duration.count() << " [s]";
      if(sample_strategy_->reestimate_hyperparameters(iteration)) {
	CollapsedGibbsDMC* which_sampler = this;
	t_usage_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_);
	slot_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_);
	gov_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_gov_);
	rel_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_rel_);
      }
      if(sample_strategy_->reestimate_gov(iteration)) {
	gov_dmc_.reestimate_collapsed_parameters(c_template_gov_);
	reestimated.g = true;
      }
      if(sample_strategy_->reestimate_rel(iteration)) {
	rel_dmc_.reestimate_collapsed_parameters(c_slot_rel_);
	reestimated.r = true;
      }
      if(sample_strategy_->reestimate_template_usage(iteration)) {
	t_usage_dmc_.reestimate_collapsed_parameters(c_template_);
	reestimated.tu = true;

	ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw_ptr =
	  dynamic_cast< ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* >
	  (
	   sit_lab()
	   );
	if(! tsw_ptr ) {
	  ERROR << "Could not get a ThriftSmartWriter; skipping labeling";
	} else {
	  // This needs to be in here for multithreaded
	  //ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw1 = tsw_ptr->clone();
	  //tsw1->reset_base(tsw1->form_id( doc.id ));
	  for(int di = 0; di < num_docs_; di++) {
	    DocT doc;
	    corpus->fill(di, &doc);
	    tsw_ptr->reset_base(tsw_ptr->form_id( doc.id ));
	    std::string suffix = "concrete::CommunicationTagging::epoch" + std::to_string(epoch) + "::iter" + std::to_string(iteration);
	    tsw_ptr->advance(suffix);
	    write_usage_posterior(tsw_ptr, t_usage_dmc_.collapsed_params(di));
	  }
	  //delete tsw1;
	}

      }
      if(sample_strategy_->reestimate_slot_usage(iteration)) {
	slot_dmc_.reestimate_collapsed_parameters(c_slot_);
	reestimated.su = true;
      }
      // print template-governors
      if(reestimated.g && reestimated.r && reestimated.tu && reestimated.su &&
	 (iteration+1) % print_struct.print.tg == 0) {
	print_dist_diag(print_struct, 10);
      }
      // // print template usage
      // if((iteration+1) % sw_wrapper.print.usage_t == 0) {
      // }
    }
  }

  template <typename CorpusT>
  void CollapsedGibbsDMC::print_labeling(CorpusT* corpus, const std::string& f_name, const VocabT& gov_vocab, const VocabT& rel_vocab) {
    std::ofstream myfile;
    myfile.open(f_name);
    myfile << "doc_id\tent\ttemplate\tslot\tmention\tgov_o\trel_o\tgov_k\trel_k\tlatent_kinds\n";
    //myfile << "doc_id\tent\ttemplate\tslot\tmention\tgov_o\trel_o\tgov_v\trel_v\tgov_k\trel_k\tlatent_kinds\n";
    for(int di = 0; di < num_docs_; di++) {
      DocT doc;
      corpus->fill(di, &doc);
      const int num_ents = doc.entities.size();
      for(int ei = 0; ei < num_ents; ei++) {
	const int prev_template = template_assignments_[di][ei];
	const int prev_slot = slot_assignments_[di][ei];
	std::string doc_ent_str = doc.id + "\t" + std::to_string(ei) + 
	  "\t" + std::to_string(prev_template) + "\t" + std::to_string(prev_slot) + "\t";
	const auto& entity = doc.entities[ei];
	const int num_ments = entity.mentions.size();
	size_t mi = 0;
	for(const auto& mention : entity.mentions) {
	  for(const auto& pa : mention.structures) {
	    if(pa.annot_level != annot_level_) continue;
	    std::string actual_gov = gov_vocab.word(pa.predicate.word);
	    std::string actual_rel = rel_vocab.word(pa.relation);
	    myfile << doc_ent_str;
	    myfile << mi << "\t" << actual_gov << "\t" << actual_rel << "\t"; 
	    //myfile << mention.gov().view() << "\t" << mention.rel() << "\t";
	    myfile << "-1\t-1\tfalse\n";
	  }
	  ++mi;
	}
      }
    }
    myfile.close();
    INFO << "Wrote labeling to " << f_name;
  }

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  template <typename CorpusT>
  void CollapsedGibbsDMCWithKinds::init(int num_templates, 
					const std::vector<int>& num_slots,
					CorpusT* corpus,
					const VocabF& gov_kind_vocab,
					const VocabF& rel_kind_vocab,
					const VocabT& gov_vocab,
					const VocabT& rel_vocab) {
    num_docs_ = corpus->num_docs();
    template_assignments_ = std::vector<std::vector<int> >(num_docs_);
    slot_assignments_ = std::vector<std::vector<int> >(num_docs_);
    num_slots_ = num_slots;
    if(latent_kinds_) {
      gov_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
      rel_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
    }
    if(num_gov_kinds_ < 0 || num_rel_kinds_ < 0) {
      ERROR << "There was an issue in getting either the number of latent governor kinds ( " << num_gov_kinds_ << ") or the number of latent relation kinds (" << num_rel_kinds_ << "). Please fix.";
      throw 10;
    }

    t_usage_dmc_ = dmc::gdmc(num_templates_,
			     std::vector<double>(num_templates_, shp_.h_theta),
			     num_docs_);
    slot_dmc_ = dmc::gdmc(num_slots_[0],
			  std::vector<double>(num_slots_[0], shp_.h_slot),
			  num_templates_);
    gov_kind_gov_dmc_ = dmc::gdmc(gov_vocab.num_words(),
				  std::vector<double>(gov_vocab.num_words(), shp_.h_gov),
				  num_gov_kinds_);
    rel_kind_rel_dmc_ = dmc::gdmc(rel_vocab.num_words(),
				  std::vector<double>(rel_vocab.num_words(), shp_.h_rel),
				  num_rel_kinds_);
    gov_kind_dmc_ = dmc::gdmc(num_gov_kinds_,
			      std::vector<double>(num_gov_kinds_, shp_.h_gov_kind),
			      num_templates_);
    rel_kind_dmc_ = dmc::mfgdmc(num_rel_kinds_,
				std::vector<double>(num_rel_kinds_, shp_.h_rel_kind),
				num_templates_,
				num_slots_);

    c_template_sum_.resize(num_templates_, 0);
    for(int gk = 0; gk < num_gov_kinds_; ++gk) {
      c_gov_kind_gov_.push_back(std::vector<int>(gov_vocab.num_words(), 0));
      c_gov_kind_gov_sum_.push_back(0);
    }
    for(int rk = 0; rk < num_rel_kinds_; ++rk) {
      c_rel_kind_rel_.push_back(std::vector<int>(rel_vocab.num_words(), 0));
      c_rel_kind_rel_sum_.push_back(0);
    }
    for(int t = 0; t < num_templates_; t++) {
      c_template_gov_kind_.push_back(std::vector<int>(num_gov_kinds_, 0));
      c_slot_.push_back(std::vector< int>(num_slots_[t]));
      c_slot_rel_kind_.push_back(std::vector< std::vector< int> >(num_slots_[t]));
      for(int s = 0; s < num_slots_[t]; s++) {
	c_slot_rel_kind_[t][s] = std::vector<int>(num_rel_kinds_);
      }
    }
    for(int di = 0; di < num_docs_; di++) {
      c_template_.push_back(std::vector<int>(num_templates_, 0));
      DocT doc;
      corpus->fill(di, &doc);
      const int num_ents = doc.entities.size();
      //template_assignments_.push_back(std::vector<int>(num_ents,0));
      //slot_assignments_.push_back(std::vector<int>(num_ents,0));
      template_assignments_[di] = std::vector<int>(num_ents, 0);
      slot_assignments_[di]     = std::vector<int>(num_ents, 0);
      if( latent_kinds_ ) {
	gov_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
	rel_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
      }
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
	auto& entity = doc.entities[ei];

	// if the kinds are NOT latent, then we just read off kind/lemma
	// stats from the mentions themselves
	// otherwise, if the kinds are latent, then we'll need to acquire 
	// stats about them in order to sample them later
	if(!latent_kinds_) {
	  // and update the possibly-observed/observed
	  for(const auto& mention : entity.mentions) {
	    int syn_level = minsky::find_level(mention, minsky::AnnotationLevel::SYNTAX);
	    int sem_level = minsky::find_level(mention, minsky::AnnotationLevel::SEMANTIC);
	    if(syn_level < 0 || sem_level < 0) {
	      ERROR << "Did not expect to find syntactic level " << syn_level << " and semantic level " << sem_level;
	      throw 10;
	    }
	    const int gk = mention.structures[sem_level].predicate.word;
	    const int gs = mention.structures[syn_level].predicate.word;
	    c_template_gov_kind_[sampled_template][gk] += 1;
	    c_gov_kind_gov_sum_[gk] += 1;
	    c_gov_kind_gov_[gk][gs] += 1;
	    // recall: we already updated the marginal c_gov_kind_gov_sum_ above

	    const int rk = mention.structures[sem_level].relation;
	    const int rs = mention.structures[syn_level].relation;
	    c_slot_rel_kind_[sampled_template][sampled_slot][rk] += 1;
	    c_rel_kind_rel_sum_[rk] += 1;
	    c_rel_kind_rel_[rk][rs] += 1;
	    // recall: we already updated the marginal c_rel_kind_rel_sum_ above
	  }
	} else {
	  const int num_ments = entity.mentions.size();
	  gov_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));
	  rel_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));	  
	  for(int mi = 0; mi < num_ments; ++mi) {
	    const auto& mention = entity.mentions[mi];
	    int syn_level = minsky::find_level(mention, minsky::AnnotationLevel::SYNTAX);
	    if(syn_level < 0) {
	      ERROR << "Did not expect to find syntactic level " << syn_level;
	      throw 10;
	    }
	    const int gv = mention.structures[syn_level].predicate.word;
	    int sampled_gov_kind = gsl_rng_uniform_int(rnd_gen_, num_gov_kinds_);
	    ++c_template_gov_kind_[sampled_template][sampled_gov_kind];
	    ++c_gov_kind_gov_sum_[sampled_gov_kind];
	    gov_kind_assignments_[di][ei][mi] = sampled_gov_kind;
	    ++c_gov_kind_gov_[sampled_gov_kind][gv];
	    ////////////////////////////
	    int sampled_rel_kind = gsl_rng_uniform_int(rnd_gen_, num_rel_kinds_);
	    ++c_slot_rel_kind_[sampled_template][sampled_slot][sampled_rel_kind];
	    ++c_rel_kind_rel_sum_[sampled_rel_kind];
	    rel_kind_assignments_[di][ei][mi] = sampled_rel_kind;
	    const int rv = mention.structures[syn_level].relation;
	    ++c_rel_kind_rel_[sampled_rel_kind][rv];
	  }
	}
      }
    }
  }

  template <typename CorpusT>
  void CollapsedGibbsDMCWithKinds::init(int num_templates, int num_slots,
					CorpusT* corpus,
					const VocabF& gov_kind_vocab,
					const VocabF& rel_kind_vocab,
					const VocabT& gov_vocab,
					const VocabT& rel_vocab) {
    init(num_templates_, std::vector<int>(num_templates_, num_slots),
	 corpus, gov_kind_vocab, rel_kind_vocab,
	 gov_vocab, rel_vocab);
  }

  template <typename CorpusT>
  void CollapsedGibbsDMCWithKinds::init(CorpusT* corpus,
					const VocabF& gov_kind_vocab,
					const VocabF& rel_kind_vocab,
					const VocabT& gov_vocab,
					const VocabT& rel_vocab) {
    init(num_templates_, num_slots_,
	 corpus, gov_kind_vocab, rel_kind_vocab,
	 gov_vocab, rel_vocab);
  }

  // template <typename G, typename R> void heldout_init(int num_templates, 
  // 						      const std::vector<int>& num_slots,
  // 						      const Vocabulary<G>& gov_kind_vocab,
  // 						      const Vocabulary<R>& rel_kind_vocab) {
  //   if(latent_kinds_) {
  //     gov_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
  //     rel_kind_assignments_ = std::vector< std::vector< std::vector<int> > >(num_docs_);
  //   }

  //   num_gov_kinds_ = model_->num_gov_kinds();
  //   num_rel_kinds_ = model_->num_rel_kinds();
  //   if(num_gov_kinds_ < 0 || num_rel_kinds_ < 0) {
  //     ERROR << "There was an issue in getting either the number of latent governor kinds ( " << num_gov_kinds_ << ") or the number of latent relation kinds (" << num_rel_kinds_ << "). Please fix.";
  //     throw 10;
  //   }

  //   t_usage_dmc_ = dmc::gdmc(model_->num_templates(), model_->hyper_theta(), num_docs_);

  //   // since this is heldout, we don't re-initialize any other DMC params
  //   c_gov_kind_gov_.clear();
  //   c_gov_kind_gov_sum_.clear();
  //   c_template_sum_.resize(num_templates, 0);
  //   for(int gk = 0; gk < gov_kind_vocab.num_words(); ++gk) {
  //     c_gov_kind_gov_.push_back(std::vector<int>(actual_gov_vocab_.num_words(), 0));
  //     c_gov_kind_gov_sum_.push_back(0);
  //   }
  //   c_rel_kind_rel_.clear();
  //   c_rel_kind_rel_sum_.clear();
  //   for(int rk = 0; rk < rel_kind_vocab.num_words(); ++rk) {
  //     c_rel_kind_rel_.push_back(std::vector<int>(actual_rel_vocab_.num_words(), 0));
  //     c_rel_kind_rel_sum_.push_back(0);
  //   }
  //   c_template_gov_kind_.clear();
  //   c_slot_.clear();
  //   c_slot_rel_kind_.clear();
  //   for(int t = 0; t < num_templates; t++) {
  //     c_template_gov_kind_.push_back(std::vector<int>(gov_kind_vocab.num_words(), 0));
  //     c_slot_.push_back(std::vector< int>(num_slots[t]));
  //     c_slot_rel_kind_.push_back(std::vector< std::vector< int> >(num_slots[t]));
  //     num_slots_[t] = num_slots[t];
  //     for(int s = 0; s < num_slots_[t]; s++) {
  // 	c_slot_rel_kind_[t][s] = std::vector<int>(rel_kind_vocab.num_words());
  //     }
  //   }
  //   c_template_.clear();
  //   template_assignments_.clear();
  //   slot_assignments_.clear();
  //   for(int di = 0; di < num_docs_; di++) {
  //     c_template_.push_back(std::vector<int>(num_templates, 0));
  //     const D& doc = (*corpus_)[di];
  //     const int num_ents = num_entities_[di];
  //     template_assignments_.push_back(std::vector<int>(num_ents, 0));
  //     slot_assignments_.push_back(std::vector<int>(num_ents, 0));
  //     //   	template_assignments_[di] = std::vector<int>(num_ents, 0);
  //     // slot_assignments_[di]     = std::vector<int>(num_ents, 0);

  //     if( latent_kinds_ ) {
  // 	gov_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
  // 	rel_kind_assignments_.push_back(std::vector< std::vector< int > >(num_ents));
  //     }
  //     for(int ei = 0; ei < num_ents; ei++){
  // 	int sampled_template = gsl_rng_uniform_int(rnd_gen_, num_templates);
  // 	template_assignments_[di][ei] = sampled_template;
  // 	int sampled_slot = gsl_rng_uniform_int(rnd_gen_, num_slots_[sampled_template]);
  // 	slot_assignments_[di][ei] = sampled_slot;
  // 	// update count tables
  // 	// these are the latent/latent count tables
  // 	++c_template_[di][sampled_template];
  // 	++c_template_sum_[sampled_template];
  // 	++c_slot_[sampled_template][sampled_slot];
  // 	// now update the latent/possibly-observed (kind) count tables
  // 	auto& entity = doc[ei];

  // 	// if the kinds are NOT latent, then we just read off kind/lemma
  // 	// stats from the mentions themselves
  // 	// otherwise, if the kinds are latent, then we'll need to acquire 
  // 	// stats about them in order to sample them later
  // 	if(!latent_kinds_) {
  // 	  std::map<G, int> gov_kind_hist = entity.gov_histogram();
  // 	  for(auto ghiter : gov_kind_hist) {
  // 	    int gi = gov_kind_vocab.index(ghiter.first);
  // 	    c_template_gov_kind_[sampled_template][gi] += ghiter.second;
  // 	    c_gov_kind_gov_sum_[gi] += ghiter.second;
  // 	  }
  // 	  std::map<G, int> rel_kind_hist = entity.rel_histogram();
  // 	  for(auto rhiter : rel_kind_hist) {
  // 	    int ri = rel_kind_vocab.index(rhiter.first);
  // 	    c_slot_rel_kind_[sampled_template][sampled_slot][ri] += rhiter.second;
  // 	    c_rel_kind_rel_sum_[ri] += rhiter.second;
  // 	  }
  // 	  // and update the possibly-observed/observed
  // 	  ferrum::pair_icount<int> gov_hist = entity.gov_view_lemma_histogram(gov_kind_vocab, actual_gov_vocab_);
  // 	  for(auto ghiter : gov_hist) {
  // 	    std::pair<int,int> p = ghiter.first;
  // 	    int count = ghiter.second;
  // 	    c_gov_kind_gov_[p.first][p.second] += count;
  // 	    // recall: we already updated the marginal c_gov_kind_gov_sum_ above
  // 	  }
  // 	  ferrum::pair_icount<int> rel_hist = entity.rel_view_str_histogram(rel_kind_vocab, actual_rel_vocab_);
  // 	  for(auto rhiter : rel_hist) {
  // 	    std::pair<int,int> p = rhiter.first;
  // 	    int count = rhiter.second;
  // 	    c_rel_kind_rel_[p.first][p.second] += count;
  // 	    // recall: we already updated the marginal c_rel_kind_rel_sum_ above
  // 	  }
  // 	} else {
  // 	  const int num_ments = entity.num_mentions();
  // 	  gov_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));
  // 	  rel_kind_assignments_[di].push_back(std::vector< int >(num_ments, 0));	  
  // 	  for(int mi = 0; mi < num_ments; ++mi) {
  // 	    const auto& mention = entity.mention(mi);
  // 	    // int gkv = gov_kind_vocab.index(mention.gov().view());
  // 	    int sampled_gov_kind = gsl_rng_uniform_int(rnd_gen_, num_gov_kinds_);
  // 	    ++c_template_gov_kind_[sampled_template][sampled_gov_kind];
  // 	    ++c_gov_kind_gov_sum_[sampled_gov_kind];
  // 	    gov_kind_assignments_[di][ei][mi] = sampled_gov_kind;
  // 	    int gv = actual_gov_vocab_.index(mention.gov().lemma());
  // 	    ++c_gov_kind_gov_[sampled_gov_kind][gv];
  // 	    ////////////////////////////
  // 	    int sampled_rel_kind = gsl_rng_uniform_int(rnd_gen_, num_rel_kinds_);
  // 	    ++c_slot_rel_kind_[sampled_template][sampled_slot][sampled_rel_kind];
  // 	    ++c_rel_kind_rel_sum_[sampled_rel_kind];
  // 	    rel_kind_assignments_[di][ei][mi] = sampled_rel_kind;
  // 	    // int rkv = rel_kind_vocab.index(mention.rel());
  // 	    int rv = actual_rel_vocab_.index(mention.rel_str());
  // 	    ++c_rel_kind_rel_[sampled_rel_kind][rv];
  // 	  }
  // 	}
  //     }
  //   }
  // }

  // FIX & MOVE
  template <typename CorpusT, template <typename> class TSW, typename... TSWArgs >
  void CollapsedGibbsDMCWithKinds::learn
  (
   CorpusT* corpus,
   int epoch,
   const StringDiscreteKindPrinter& print_struct,
   DKVWriters* sw_wrapper,
   bool heldout,
   TSWArgs... tsw_args
   ) {
    ferrum::Timer timer("Collapsed Gibbs");
    reestimated reest = {heldout, heldout, heldout, heldout, false, heldout};
    bool do_labeling = print_struct.print.usage_t > 0;
    ferrum::ConcreteSituationLabeler<ThriftProtocol, TSW> sit_lab(do_labeling,
								  print_struct.print.usage_t);
    sit_lab.template make_with_args< std::string, TSWArgs... >("template_props.tcompact", tsw_args...);
    sit_lab.corpus(corpus);
    for(int iteration = 0; iteration < sample_strategy_->num_iterations; ++iteration) {
      for(int di = 0; di < num_docs_; di++) {
	DocT doc;
	corpus->fill(di, &doc);
	const int num_ents = doc.entities.size();
	for(int ei = 0; ei < num_ents; ei++) {
	  bool floating_hists_updated = false;
	  const auto& entity = doc.entities[ei];
	  if(sample_strategy_->sample_template(iteration, di, ei)) {
	    update_floating_hists(entity, di, ei);
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
	      update_floating_hists(entity, di, ei);
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
	  if(latent_kinds_ || heldout) {
	    const auto& entity = doc.entities[ei];
	    const int num_ments = entity.mentions.size();
	    for(int mi = 0; mi < num_ments; ++mi) {
	      const auto& mention = entity.mentions[mi];
	      int which = minsky::find_level(mention, minsky::AnnotationLevel::SEMANTIC);
	      if(which < 0) {
		ERROR << "Did not expect to find -1 level; throwing";
		throw 10;
	      }
	      floating_mention_obs_.first =  mention.structures[which].predicate.word;
	      floating_mention_obs_.second = mention.structures[which].relation;
	      if(sample_strategy_->sample_gov_kind(iteration, di, ei, mi, latent_kinds_)) {
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
	      if(sample_strategy_->sample_rel_kind(iteration, di, ei, mi, latent_kinds_)) {
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
      }
      auto duration = timer.lap();
      INFO << "Iteration " << (iteration+1) << " took " << duration.count() << " [s]";
      if(sample_strategy_->reestimate_hyperparameters(iteration)) {
	CollapsedGibbsDMCWithKinds* which_sampler =  this;
	  //heldout_sampler != NULL ? heldout_sampler : this;
	// if(heldout_sampler != NULL) {
	//   heldout_sampler->learn(gov_kind_vocab, rel_kind_vocab);
	// }
	t_usage_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_);
	//slot_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_);
	gov_kind_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_template_gov_kind_);
	rel_kind_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_slot_rel_kind_);
	gov_kind_gov_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_gov_kind_gov_);
	rel_kind_rel_dmc_.reestimate_hyperparameters_wallach1(which_sampler->c_rel_kind_rel_);
      }
      attempt_reestimation(iteration, reest);
      if(sample_strategy_->reestimate_template_usage(iteration)) {
	ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw_ptr =
	  dynamic_cast< ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* >
	  (
	   sit_lab()
	   );
	if(! tsw_ptr ) {
	  ERROR << "Could not get a ThriftSmartWriter; skipping labeling";
	} else {
	  for(int di = 0; di < num_docs_; di++) {
	    DocT doc;
	    corpus->fill(di, &doc);
	    tsw_ptr->reset_base(tsw_ptr->form_id( doc.id ));
	    std::string suffix = "concrete::CommunicationTagging::epoch" + std::to_string(epoch) + "::iter" + std::to_string(iteration);
	    tsw_ptr->advance(suffix);
	    write_usage_posterior(tsw_ptr, t_usage_dmc_.collapsed_params(di));
	  }
	}
      }

      if(reest.tu &&
	 (heldout || (reest.g && reest.r &&
		      reest.su && reest.gk &&
		      reest.rk)) &&
	 (iteration+1) % print_struct.print.tg == 0) {
	print_dist_diag(print_struct, 10);
      }
    } // end for(iteration)
  } // end learn method
  
  template <typename CorpusT>
  void CollapsedGibbsDMCWithKinds::print_labeling
  (CorpusT* corpus, const std::string& f_name,
   const VocabF& gkv, const VocabF& rkv,
   const VocabT& gv, const VocabT& rv) {
    std::ofstream myfile;
    myfile.open(f_name);
    myfile << "doc_id\tent\ttemplate\tslot\tmention\tgov_o\trel_o\tgov_k\trel_k\tlatent_kinds\n";
    for(int di = 0; di < num_docs_; di++) {
      DocT doc;
      corpus->fill(di, &doc);
      const int num_ents = doc.entities.size();
      for(int ei = 0; ei < num_ents; ei++) {
	const int a_template = template_assignments_[di][ei];
	const int a_slot = slot_assignments_[di][ei];
	std::string doc_ent_str = doc.id + "\t" + std::to_string(ei) + 
	  "\t" + std::to_string(a_template) + "\t" + std::to_string(a_slot) + "\t";
	const auto& entity = doc.entities[ei];
	const int num_ments = entity.mentions.size();
	size_t mi = 0;
	for(const auto& mention : entity.mentions) {
	  std::string gov, rel, govk, relk;
	  for(const auto& pa : mention.structures) {
	    if(pa.annot_level == minsky::AnnotationLevel::SYNTAX) {
	      gov = gv.word(pa.predicate.word);
	      rel = rv.word(pa.relation);
	    } else if(pa.annot_level == minsky::AnnotationLevel::SEMANTIC) {
	      if(latent_kinds_) {
		govk = std::to_string(gov_kind_assignments_[di][ei][mi]);
		relk = std::to_string(rel_kind_assignments_[di][ei][mi]);
	      } else {
		govk = gkv.word(pa.predicate.word);
		relk = rkv.word(pa.relation);
	      }
	    }
	  }
	  myfile << doc_ent_str << mi << "\t";
	  myfile << gov << "\t" << rel << "\t";
	  myfile << govk << "\t" << relk << "\t";
	  myfile << (latent_kinds_ ? "true" : "false") << "\n";
	  ++mi;
	} // end mention loop
      } // end entity loop
    }
    myfile.close();
    INFO << "Wrote labeling to " << f_name;
  }

} // end namespace ferrum

#endif // FERRUM_LIBNAR_CRTLDA_SAMPLING_TCC_
