#include "ferrum/sage_defs.hpp"
#include "ferrum/sage_tm_variational.hpp"
#include "ferrum/redis_corpus.hpp"
#include "ferrum/stats.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/timer.hpp"
#include "ferrum/version.hpp"

#include <chrono>
#include <utility>

namespace ferrum {
  SageTMVariational::SageTMVariational() :
    num_docs_(0), num_topics_(0), initialized_(false) {
  }
  SageTMVariational::SageTMVariational
  (
   int nt,
   const std::vector<double>& usage_h,
   const std::vector<double>& word_h
   ) :
    num_docs_(0),
    num_topics_(nt),
    topic_usage_hypers_(usage_h),
    word_hypers_(word_h),
    initialized_(false) {
  }
  SageTMVariational::SageTMVariational
  (
   int nt,
   int vocab_size,
   double usage_h,
   double word_h
   ) :
    num_docs_(0),
    num_topics_(nt),
    topic_usage_hypers_( Vector1D(num_topics_, usage_h) ),
    word_hypers_( Vector1D(vocab_size, word_h) ),
    initialized_(false) {
  }
  SageTMVariational::~SageTMVariational() {
  }

  void SageTMVariational::num_topics(int nt) {
    num_topics_ = nt;
  }
  int SageTMVariational::num_topics() {
    return num_topics_;
  }

  // TODO: probably refactor this out w/ the other *sage* inference methods
  typename SageTMVariational::Vector2D
  SageTMVariational::get_distributions(const Vector2D* ptr) {
    Vector2D topics;
    typedef Vector2D::value_type Inner;
    for(const auto& utopic : *ptr) {
      const double norm = ferrum::sum(utopic);
      Inner topic(ferrum::scalar_product(1.0/norm, utopic));
      topics.push_back(topic);
    }
    return topics;
  }
  typename SageTMVariational::Vector2D
  SageTMVariational::get_distributions(std::vector< typename SageTMVariational::STopicType >* vec) {
    typedef std::vector< double > VecType;
    std::vector< VecType > rvec;
    for(const STopicType& topic : *vec) {
      VecType sub = topic.template eta_as<VecType>(true);
      rvec.push_back( std::move( sub ) );
    }
    return rvec;
  }

  void SageTMVariational::update_topic_assignments
  (
   Vector1D* assign_params,
   const Vector1D& usage_params,
   const Vector1D& expected_word_counts
   ) {
    ferrum::copy(&usage_params, assign_params);
    ferrum::sum_in_first(assign_params, expected_word_counts);
    ferrum::prob_from_unnorm_lp(assign_params);
  }

  void SageTMVariational::update_topic_usage(std::vector<double>* var_assign,
					     const std::vector<double>& buffer) {
    *var_assign = topic_usage_hypers_;
    ferrum::sum_in_first(var_assign, buffer);
  }

  // SAGE-specific
  void SageTMVariational::_update_topic_word_usage(int index_topic, double i1, double i2) {
    INFO << "Beginning to fit topic/word distribution " << index_topic;
    STopicType& sage_topic = var_topic_word_params_[index_topic];
    std::vector<double>& usable_counts = buffer_topic_word_params_[index_topic];
    int opt_status = sage_topic.fit_topic_svi(&usable_counts, 1.0, i1, i2);
    INFO << "Topic " << index_topic << " has optimization status " << opt_status;
    // std::vector<DoublePair>& dp = var_variance_template_verb_params_[index_template];
    // sage_topic.update_variances(&dp, i1, i2);
  }
  void SageTMVariational::update_topic_word_usage(int num_threads, double i1, double i2) {
    omp_set_num_threads(num_threads);
#pragma omp parallel for
    for(int ti = 0; ti < num_topics_; ++ti) {
      _update_topic_word_usage(ti, i1, i2);
    }
  }

  // template <typename P>
  // void SageTMVariational::write_situations
  // (
  //  const DType& doc,
  //  ferrum::thrift::ThriftSmartWriter<P>* tsw,
  //  const Vector2D& t_assign,
  //  const Vector2D& s_assign
  //  ) {
  //   assert(ferrum::num_entities(doc) == t_assign.size());
  //   assert(t_assign.size() == s_assign.size());
  //   // this is currently a "mode" serializer
  //   Vector1D t_max = ferrum::row_arg_max(t_assign);
  //   Vector1D s_max = ferrum::row_arg_max(s_assign);
  //   concrete::util::uuid_factory uf;
  //   size_t i = 0;
  //   for(const auto& entity : get_entities(doc) ) {
  //     concrete::Situation situation;
  //     // required: Situation.uuid
  //     {
  // 	concrete::UUID uuid;
  // 	uuid.__set_uuidString(uf.get_uuid());
  // 	situation.__set_uuid(uuid);
  //     }
  //     // required: Situation.situationType
  //     {
  // 	situation.__set_situationType("EVENT");
  //     }
  //     std::string situation_kind = std::to_string(t_max[i]);
  //     // optional: Situation.situationKind
  //     {
  // 	situation.__set_situationKind(situation_kind);
  //     }
  //     std::string situation_role = std::to_string(s_max[i]);
  //     // now iterate over the mentions
  //     std::vector<concrete::UUID> si_mention_id_list;
  //     for(const auto& doc_mention : get_mentions(entity) ) {
  // 	//const auto& conc_mention;
  // 	concrete::SituationMention sm;
  // 	// required: uuid
  // 	{
  // 	  concrete::UUID uuid;
  // 	  uuid.__set_uuidString(uf.get_uuid());
  // 	  sm.__set_uuid(uuid);
  // 	  si_mention_id_list.push_back(uuid);
  // 	}
  // 	// optional: situationKind == template assignment
  // 	{
  // 	  sm.__set_situationKind(situation_kind);
  // 	}
  // 	// fill MentionArgument
  // 	{
  // 	  concrete::MentionArgument ma;
  // 	  // optional: MentionArgument.role == slot assignment
  // 	  ma.__set_role(situation_role);
  // 	  concrete::UUID ma_ent_uuid;
  // 	  if(doc_mention.id.size() == 0) {
  // 	    WARN << "Document [" << doc.id << "] has entity [" << entity.id << "] with an unset mention.id";
  // 	  }
  // 	  ma_ent_uuid.__set_uuidString(doc_mention.id);
  // 	  ma.__set_entityMentionId(ma_ent_uuid);
  // 	  // it would be nice to set the .tokens field too... oh well
  // 	  // required: argumentList
  // 	  {
  // 	    std::vector<concrete::MentionArgument> sm_args;
  // 	    sm_args.push_back(ma);
  // 	    sm.__set_argumentList(sm_args);
  // 	  }

  // 	  // fill Situation.Argument
  // 	  concrete::Argument sit_arg;
  // 	  sit_arg.__set_role(situation_role);
  // 	  sit_arg.__set_entityId(ma_ent_uuid);
  // 	}
  // 	// optional: situationType
  // 	{
  // 	  sm.__set_situationType("EVENT");
  // 	}
  // 	ferrum::thrift::save<P, concrete::SituationMention>(tsw, sm);
  //     }
  //     if(si_mention_id_list.size() > 0) {
  // 	situation.__set_mentionIdList( si_mention_id_list );
  //     }
  //     ferrum::thrift::save<P, concrete::Situation>(tsw, situation);
  //     ++i;
  //   }
  // }
  template <typename P>
  void SageTMVariational::write_usage_posterior
  (
   ferrum::thrift::ThriftSmartWriter<P>* tsw,
   const Vector1D& usage_post
   ) {
    concrete::CommunicationTagging ct;
    concrete::util::uuid_factory uf;
    // required: uuid
    concrete::util::add_uuid(ct, uf);
    // required: metadata
    concrete::util::add_metadata(ct,
				 "SAGE topic induction: " + ferrum::FERRUM_GIT_SHA);
    ct.__set_taggingType("sage-topics");
    // loop through a (sorted) posterior, setting both the tagList and confidenceList
    const double norm = ferrum::sum(usage_post);
    std::vector<size_t> sorted_topic = ferrum::sort_indices(usage_post, false);
    ct.__isset.tagList = true;
    ct.__isset.confidenceList = true;
    ct.tagList.resize(num_topics_);
    ct.confidenceList.resize(num_topics_);
    std::vector<std::string>& tl = ct.tagList;
    std::vector<double>& cl = ct.confidenceList;
    for(size_t i = 0; i < (size_t)num_topics_; ++i) {
      size_t which = sorted_topic[i];
      tl[i] = std::to_string(which);
      cl[i] = (double)usage_post[which] / norm ;
    }
    ferrum::thrift::save<P, concrete::CommunicationTagging>(tsw, ct);
  }

  template <typename CorpusT>
  double SageTMVariational::e_step
  (
   const VStrategy& strategy,
   int learn_iter,
   const size_t batch_start,
   const size_t batch_end,
   BaseSituationLabeler* sit_lab,
   const size_t offset,
   const int epoch
   ) {
    ferrum::Timer estep_timer(__func__);
    // TODO: there's *GOT* to be a better way of handling these aggregated values
    // but for now, we'll just eat the memory cost (which will be pretty small)
    Vector2D words_grads(num_words_, Vector1D(num_topics_, 0.0));
    if(strategy.em_verbosity >= 1) {
      INFO << "E-step";
    }
    renormalize_sage_distributions();
    //double elbo_base = 0.0; // cache parts of the elbo computation
    for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
      ferrum::set(&(buffer_topic_word_params_[ti]), 0.0);
    }
    for(size_t vi = 0; vi < (size_t)num_words_; ++vi) {
      for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
	words_grads[vi][ti] = var_topic_word_params_[ti].l_probability(vi);
	//elbo_base -= stats::gamma::entropy(var_variance_template_verb_params_[ti][vi]);
      }
    }
    const size_t batch_size = batch_end - batch_start;
    omp_set_num_threads(strategy.num_e_threads);
    double elbo = 0.0;
#pragma omp parallel for \
  reduction(+:elbo)
    for(size_t di = 0; di < batch_size; ++di) {
      //double doc_elbo = 0.0;//, prev_elbo = std::numeric_limits<double>::max();
      const size_t num_words = words_in_docs_[di].size();
      Vector2D var_topic_assign_params(num_words);
      // These two matrices contain the accumulated gradient-log-normalizers
      Vector2D expected_word_counts(num_words, Vector1D(num_topics_, 0.0));
      // thread-specific buffer for doc-topic usage
      Vector1D butp(num_topics_, 0.0);
      CorpusT* corp_for_labeling = NULL;
      typename CorpusT::DocType doc;
      if(sit_lab != NULL && sit_lab->perform_labeling()) {
	corp_for_labeling = sit_lab->corpus< CorpusT* >();
	assert(corp_for_labeling != NULL);
	//doc = corp_for_labeling->operator[](di + offset);
	corp_for_labeling->fill(di+offset, &doc);
      }
      for(size_t word_idx = 0; word_idx < num_words; ++word_idx) {
	var_topic_assign_params[word_idx] = vi_init_.assignment(num_topics_);
	// int word_count = word_type_counts_[di][word_idx];
	// int obs_word = words_in_docs_[di][word_idx];
	// auto fg_col = ferrum::column(words_grads, obs_word);
	// ferrum::linear_combination_in_first
	//   (
	//    &(expected_word_counts[word_idx]),
	//    fg_col,
	//    1.0,
	//    (double)word_count
	//    );
      }
      double elbo_change = 0.0;
      Vector2D v_select_elbo;
      const bool compute_elbo = (strategy.compute_elbo || strategy.elbo_as_break);
      // if(compute_elbo) {
      // 	for(size_t ent_idx = 0; ent_idx < num_ents; ++ent_idx) {
      // 	  v_select_elbo.push_back(Vector1D(num_templates_));
      // 	  for(size_t ti = 0; ti < num_templates_; ++ti) {
      // 	    v_select_elbo[ent_idx][ti] =
      // 	      var_template_verb_params_[ti].elbo
      // 	      (verb_type_counts_[di][ent_idx],
      // 	       verbs_in_docs_[di][ent_idx],
      // 	       false);
      // 	  }
      // 	}
      // }
      for(int e_iter = 0; e_iter < strategy.num_e_iters; ++e_iter) {
	//	doc_elbo = elbo_base;
	if((strategy.em_verbosity >= 3) || 
	   (strategy.em_verbosity >= 2 && di % 1000 == 0)) {
	  INFO << "\tDocument " << di << " of iteration " << e_iter;
	}
	bool last_iter = 
	  ( (e_iter + 1) == strategy.num_e_iters) ||
	  ( compute_elbo && (std::abs(elbo_change) < 1E-3) );
	if(strategy.em_verbosity >= 4) {
	  INFO << "\t\tE-step sub-iteration number " << e_iter;
	}
	std::vector<double> usage_grad = 
	  dmc::dirichlet::grad_log_partition_static(var_topic_usage_params_[di]);
	for(size_t word_idx = 0; word_idx < num_words; ++word_idx) {
	  int obs_word = words_in_docs_[di][word_idx];
	  update_topic_assignments
	    (
	     &(var_topic_assign_params[word_idx]), 
	     usage_grad,
	     words_grads[obs_word] //expected_word_counts[word_idx]
	     );
	  const Vector1D& vtap = var_topic_assign_params[word_idx];
	  // const double t_cat_entropy = dmc::cat::entropy(vtap);
	  // if(compute_elbo) {
	  //   const double t_select_elbo = ferrum::dot(usage_grad, vtap);
	  //   doc_elbo += t_select_elbo - t_cat_entropy;
	  //   // compute the verb entropy, but only the "selection" elbo (no entropy)
	  //   double v_select_e = ferrum::dot(v_select_elbo[word_idx], vtap);
	  //   doc_elbo += v_select_e;
	  // }
	  ferrum::sum_in_first(&butp, vtap);
	  // and accumulate into expected counts
	  // update buffers if last time
	  if(last_iter) {
	    // update the buffers
	    {
#pragma omp critical
	      {
		for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
		  const double t_weight = vtap[ti];
		  int num_times = word_type_counts_[di][word_idx];
		  int obs = words_in_docs_[di][word_idx];
		  buffer_topic_word_params_[ti][obs] += num_times * t_weight;
		}
	      }
	    }

	    // if(sit_lab != NULL && sit_lab->perform_labeling()) {
	    //   ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw_ptr =
	    // 	dynamic_cast< ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* >
	    // 	(
	    // 	 sit_lab->operator()()
	    // 	 );
	    //   if(! tsw_ptr ) {
	    // 	ERROR << "Could not get a ThriftSmartWriter; skipping labeling";
	    //   } else {
	    // 	ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw1 = tsw_ptr->clone();
	    // 	tsw1->reset_base(tsw1->form_id( doc.id ));
	    // 	std::string suffix = "concrete::situation_labels::epoch " + std::to_string(epoch) + "::iter" + std::to_string(e_iter);
	    // 	tsw1->advance(suffix);
	    // 	write_situations(doc, tsw1, var_template_assign_params, var_slot_assign_params);
	    // 	delete tsw1;
	    //   }
	    // }
	  }
	} // end for loop(word_idx)
	update_topic_usage(&(var_topic_usage_params_[di]), butp);
	// if(compute_elbo) {
	//   doc_elbo -= dmc::dirichlet::entropy(var_template_usage_params_[di]);
	//   if(last_iter) {
	//     elbo += doc_elbo;
	//     INFO << "Doc elbo for " << di << " is " << doc_elbo;
	//   }
	// }
	ferrum::set(&butp, 0.0);
	if(sit_lab != NULL && sit_lab->perform_labeling()) {
	  ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw_ptr =
	    dynamic_cast< ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* >
	    (
	     sit_lab->operator()()
	     );
	  if(! tsw_ptr ) {
	    ERROR << "Could not get a ThriftSmartWriter; skipping labeling";
	  } else {
	    ferrum::thrift::ThriftSmartWriter<ThriftProtocol>* tsw1 = tsw_ptr->clone();
	    tsw1->reset_base(tsw1->form_id( doc.id ));
	    std::string suffix = "concrete::CommunicationTagging::epoch" + std::to_string(epoch) + "::iter" + std::to_string(e_iter);
	    tsw1->advance(suffix);
	    write_usage_posterior(tsw1, var_topic_usage_params_[di]);
	    delete tsw1;
	  }
	} // end template usage posterior labeling
      } // end for loop(e_iter)
    } // end for loop(doc \in corpus)
    return elbo;
  } // end e_step()

  const SageTMVariational::Vector1D&
  SageTMVariational::buffer_topic_word_params(size_t i) {
    return buffer_topic_word_params_[i];
  }

  void SageTMVariational::m_step(VStrategy& strategy, unsigned int batch_size) {
    ferrum::Timer mstep_timer(__func__);
    if(strategy.em_verbosity >= 1) {
      INFO << "M-step";
    }
    double rho = strategy.ssu();
    double divider = rho / (double)batch_size;
    INFO << "In m-step, interpolating with weights w1 = " << (1.0-rho) << " and w2 = " << divider;
    update_topic_word_usage(strategy.num_m_threads, 1.0 - rho, divider);
    ++(strategy.ssu);
  }

  // TODO: Fix this
  void SageTMVariational::update_hypers() {
    ERROR << "update_hypers is currently a no-op";
    // typedef dmc::DirichletVariationalClosure ClosureType;
    // ClosureType closure;
    // closure.variational_params = &var_topic_usage_params_;
    // // INFO << "The closure suff. stats are: ";
    // // ferrum::print_2d(var_topic_usage_params_);
    // double f_init = dmc::dirichlet::hyperparameters_variational_objective(usage_hypers_, &closure);
    // std::vector<double> point = dmc::dirichlet::hyperparameters_variational_nr(usage_hypers_, &closure);
    // double f_end = dmc::dirichlet::hyperparameters_variational_objective(point, &closure);
    // const double dist = ferrum::dist(usage_hypers_, point);
    // INFO << "Hyperparameter optimization moved the alpha point " << dist << " units away";
    // INFO << "Hyperparameter optimization moved the value " << (f_end - f_init) << ", from " << f_init << " to " << f_end;
    // usage_hypers_ = point;
  }

  // void SageTMVariational::_print_in_learn
  // (
  //  const VStrategy& strategy, 
  //  const StringDiscreteKindPrinter* const print_struct,
  //  DKVWriters* sw_wrapper,
  //  int epoch,
  //  int learn_iter,
  //  bool last_iter,
  //  bool model_changed
  //  ) {
  //   // TERMINAL PRINTING:: This is different than printing to file
  //   print_dist_diag(print_struct, strategy, learn_iter, last_iter);
  //   if(sw_wrapper == NULL) {
  //     return;
  //   }
  //   std::string suff = "epoch" + std::to_string(epoch) + ".emiter" + std::to_string(learn_iter);
  //   if(learn_iter % strategy.print_topics_every == 0 || last_iter) {
  //     // in order to print any of the distributions, we need the vocabs
  //     if(print_struct != NULL) {
  // 	// TODO: VOCAB DISTRIBUTION PRINTING
  // 	// TODO: make this bit more general (macro, perhaps?)
  // 	if(sw_wrapper->to_file_verb_usage()) {
  // 	  ferrum::SmartWriter* swptr = sw_wrapper->sw_verb_usage();
  // 	  std::ostream& out_stream = swptr->get(suff);
  // 	  INFO << "Writing verb distributions to " << swptr->name();
  // 	  // if we haven't updated, then we should get some form of the current topics
  // 	  auto* vocab_ptr = print_struct->g_vocab;
  // 	  if(learn_iter % strategy.update_model_every) {
  // 	    std::vector< ModelTopicType > topics;
  // 	    for(const auto& utopic : var_template_verb_params_) {
  // 	      Vector1D vec_utopic = utopic.eta_as<Vector1D>(false);
  // 	      const double norm = ferrum::sum( vec_utopic );
  // 	      ModelTopicType topic(ferrum::scalar_product(1.0/norm, vec_utopic));
  // 	      topics.push_back(topic);
  // 	    }
  // 	    model_->print_verbs(num_verbs_, *vocab_ptr, out_stream, topics);
  // 	  } else {
  // 	    model_->print_verbs(num_verbs_, *vocab_ptr, out_stream);
  // 	  }
  // 	}
  // 	if(sw_wrapper->to_model_tsv()) {
  // 	  ferrum::SmartWriter* swptr = sw_wrapper->sw_model_tsv();
  // 	  INFO << "Writing entire model as TSV to " << swptr->name();
  // 	  std::ostream& out_stream = swptr->get(suff);
  // 	  if(! (learn_iter % strategy.update_model_every)) {
  // 	    update_model(strategy.heldout);
  // 	  }
  // 	  auto cp_ps = *print_struct;
  // 	  cp_ps.num_per_gov = num_verbs_;
  // 	  cp_ps.num_per_rel = num_arcs_;
  // 	  model_->print_templates(out_stream, &cp_ps);
  // 	}
  //     } else {
  // 	WARN << "The SmartWriter wrapper was not null, but the vocab wrapper was. We cannot print any vocab distributions (akin to topics) with out.";
  //     }
  //   }
  //   // USAGE PRINTING
  //   if((learn_iter % strategy.print_usage_every == 0 && learn_iter > 0) || model_changed || last_iter) {
  //     if(sw_wrapper->to_file_template_usage()) {
  // 	ferrum::SmartWriter* swptr = sw_wrapper->sw_template_usage();
  // 	std::ostream& out_stream = swptr->get(suff);
  // 	INFO << "Writing document-template usage output to " << swptr->name();
  // 	//model_->print_template_usage(out_stream);
  // 	ferrum::print_2d_distribution(get_distributions(&var_template_usage_params_), out_stream);
  //     }
  //     if(sw_wrapper->to_file_slot_usage() ) {
  // 	ferrum::SmartWriter* swptr = sw_wrapper->sw_slot_usage();
  // 	std::ostream& out_stream = swptr->get(suff);
  // 	INFO << "Writing document-template usage output to " << swptr->name();
  // 	model_->print_slot_usage(out_stream);
  //     }
  //   }
  //   //#endif
  // }

  void SageTMVariational::renormalize_sage_distributions() {
    for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
      var_topic_word_params_[ti].renormalize();
    }
  }

  minsky::residual::ResidualTopicModel
  SageTMVariational::create_minsky_view
  (
   const ferrum::Vocabulary<std::string>& vocab
   ) {
    using namespace minsky::residual;
    using namespace minsky;
    ResidualTopicModel rtm;
    minsky::Vocab mgv = vocab.minskify();
    rtm.vocabularies.push_back(mgv);
    const int vocab_idx = 0;
    rtm.__set_background(*background_);
    for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
      Topic minsky_t;
      Frame ri = var_topic_word_params_[ti].create_minsky_frame(vocab_idx);
      minsky_t.__set_frame(ri);
      rtm.topics.push_back(minsky_t);
    }
    return rtm;
  }
} // end namespace crtlda

template
double ferrum::SageTMVariational::e_step<ferrum::InMemoryCorpus< typename ferrum::SageTMVariational::DType > >
(
 const ferrum::VStrategy& strategy,
 int learn_iter,
 const size_t batch_start,
 const size_t batch_end,
 BaseSituationLabeler* sit_lab,
 const size_t offset,
 int epoch
 );

template
double ferrum::SageTMVariational::e_step<ferrum::RedisCorpus< typename ferrum::SageTMVariational::DType > >
(
 const ferrum::VStrategy& strategy,
 int learn_iter,
 const size_t batch_start,
 const size_t batch_end,
 BaseSituationLabeler* sit_lab,
 const size_t offset,
 int epoch
 );

template
void ferrum::SageTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TBinaryProtocol>* tsw,
 const Vector1D& usage
 );
template
void ferrum::SageTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TCompactProtocol>* tsw,
 const Vector1D& usage
 );
template
void ferrum::SageTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TJSONProtocol>* tsw,
 const Vector1D& usage
 );
