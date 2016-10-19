//#include "ferrum/sage_defs.hpp"
#include "ferrum/dmc_tm_variational.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/redis_corpus.hpp"
#include "ferrum/stats.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/timer.hpp"
#include "ferrum/version.hpp"

#include <chrono>
#include <utility>

namespace ferrum {
  DMCTMVariational::DMCTMVariational() :
    num_docs_(0), num_topics_(0), initialized_(false) {
  }
  DMCTMVariational::DMCTMVariational
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
  DMCTMVariational::DMCTMVariational
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
  DMCTMVariational::~DMCTMVariational() {
  }

  void DMCTMVariational::num_topics(int nt) {
    num_topics_ = nt;
  }
  int DMCTMVariational::num_topics() {
    return num_topics_;
  }

  // // TODO: probably refactor this out w/ the other *sage* inference methods
  // typename DMCTMVariational::Vector2D
  // DMCTMVariational::get_distributions(const Vector2D* ptr) {
  //   Vector2D topics;
  //   typedef Vector2D::value_type Inner;
  //   for(const auto& utopic : *ptr) {
  //     const double norm = ferrum::sum(utopic);
  //     Inner topic(ferrum::scalar_product(1.0/norm, utopic));
  //     topics.push_back(topic);
  //   }
  //   return topics;
  // }
  // typename DMCTMVariational::Vector2D
  // DMCTMVariational::get_distributions(std::vector< typename DMCTMVariational::STopicType >* vec) {
  //   typedef std::vector< double > VecType;
  //   std::vector< VecType > rvec;
  //   for(const STopicType& topic : *vec) {
  //     VecType sub = topic.template eta_as<VecType>(true);
  //     rvec.push_back( std::move( sub ) );
  //   }
  //   return rvec;
  // }

  void DMCTMVariational::update_topic_assignments
  (
   Vector1D* assign_params,
   const Vector1D& usage_params,
   const Vector1D& expected_word_counts
   ) {
    ferrum::copy(&usage_params, assign_params);
    ferrum::sum_in_first(assign_params, expected_word_counts);
    ferrum::prob_from_unnorm_lp(assign_params);
  }

  void DMCTMVariational::update_topic_usage(std::vector<double>* var_assign,
					     const std::vector<double>& buffer) {
    *var_assign = topic_usage_hypers_;
    ferrum::sum_in_first(var_assign, buffer);
  }

  void DMCTMVariational::_update_topic_word_usage(int index_topic, double i1, double i2) {
    //INFO << "Beginning to fit topic/word distribution " << index_topic;
    STopicType& topic = var_topic_word_params_[index_topic];
    std::vector<double>& usable_counts = buffer_topic_word_params_[index_topic];
    update_global_params(word_hypers_,
			 &topic,
			 &usable_counts,
			 i1,
			 i2);
  }
  void DMCTMVariational::update_topic_word_usage(int num_threads, double i1, double i2) {
    omp_set_num_threads(num_threads);
#pragma omp parallel for
    for(int ti = 0; ti < num_topics_; ++ti) {
      _update_topic_word_usage(ti, i1, i2);
    }
  }

  template <typename P>
  void DMCTMVariational::write_usage_posterior
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
				 "DMC topic induction: " + ferrum::FERRUM_GIT_SHA);
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
  double DMCTMVariational::e_step
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
    Vector2D words_grads(num_topics_, Vector1D(num_words_, 0.0));
    if(strategy.em_verbosity >= 1) {
      INFO << "E-step";
    }
    double elbo_base = 0.0; // cache parts of the elbo computation
    (void)elbo_base;
    for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
      ferrum::set(&(buffer_topic_word_params_[ti]), 0.0);
      dmc::dirichlet::grad_log_partition_static(var_topic_word_params_[ti],
						&(words_grads[ti]));
    }
    const size_t batch_size = batch_end - batch_start;
    omp_set_num_threads(strategy.num_e_threads);
    double elbo = 0.0;
#pragma omp parallel for \
  reduction(+:elbo)
    for(size_t di = 0; di < batch_size; ++di) {
      double doc_elbo = 0.0;//, prev_elbo = std::numeric_limits<double>::max();
      const size_t num_words = words_in_docs_[di].size();
      Vector2D var_topic_assign_params(num_words);
      // These two matrices contain the accumulated gradient-log-normalizers
      //Vector2D expected_word_counts(num_words, Vector1D(num_topics_, 0.0));
      // thread-specific buffer for doc-topic usage
      Vector1D butp(num_topics_, 0.0);
      CorpusT* corp_for_labeling = corpus_for_labeling<CorpusT>(sit_lab);
      typename CorpusT::DocType doc;
      if(corp_for_labeling != NULL) {
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
      // 	for(size_t word_idx = 0; word_idx < num_words; ++word_idx) {
      // 	  v_select_elbo.push_back(Vector1D(num_topics_));
      // 	  for(size_t ti = 0; ti < num_topics_; ++ti) {
      // 	    v_select_elbo[word_idx][ti] =
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
	  auto fg_col = ferrum::column(words_grads, obs_word);
	  update_topic_assignments
	    (
	     &(var_topic_assign_params[word_idx]), 
	     usage_grad,
	     fg_col // word_grads[*,word_idx]
	     );
	  const Vector1D& vtap = var_topic_assign_params[word_idx];
	  const double t_cat_entropy = dmc::cat::entropy(vtap);
	  if(compute_elbo) {
	    const double t_select_elbo = ferrum::dot(usage_grad, vtap);
	    doc_elbo += t_select_elbo - t_cat_entropy;
	    // compute the verb entropy, but only the "selection" elbo (no entropy)
	    double v_select_e = ferrum::dot(v_select_elbo[word_idx], vtap);
	    doc_elbo += v_select_e;
	  }
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

  const DMCTMVariational::Vector1D&
  DMCTMVariational::buffer_topic_word_params(size_t i) {
    return buffer_topic_word_params_[i];
  }

  void DMCTMVariational::m_step(VStrategy& strategy, unsigned int batch_size) {
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
  void DMCTMVariational::update_hypers() {
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

  // void DMCTMVariational::_print_in_learn
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

  minsky::residual::ResidualTopicModel
  DMCTMVariational::create_minsky_view
  (
   const ferrum::Vocabulary<std::string>& vocab
   ) {
    using namespace minsky::residual;
    using namespace minsky;
    ResidualTopicModel rtm;
    minsky::Vocab mgv = vocab.minskify();
    rtm.vocabularies.push_back(mgv);
    for(size_t ti = 0; ti < (size_t)num_topics_; ++ti) {
      Topic minsky_t;
      minsky::Frame f;
      minsky::Distribution d;
      d.__set_support_size(num_words_);
      minsky::Weights w;
      STopicType cpy(var_topic_word_params_[ti]);
      ferrum::log(&cpy);
      std::copy(cpy.begin(), cpy.end(),
		std::back_inserter(w.residual));
      w.__isset.residual = true;
      d.__set_weights(w);
      d.__set_vocab_idx(0);
      f.__set_distr(d);
      minsky_t.__set_frame(f);
      rtm.topics.push_back(minsky_t);
    }
    return rtm;
  }
} // end namespace crtlda

template
double ferrum::DMCTMVariational::e_step<ferrum::RedisCorpus< typename ferrum::DMCTMVariational::DType > >
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
void ferrum::DMCTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TBinaryProtocol>* tsw,
 const Vector1D& usage
 );
template
void ferrum::DMCTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TCompactProtocol>* tsw,
 const Vector1D& usage
 );
template
void ferrum::DMCTMVariational::write_usage_posterior
(
 ferrum::thrift::ThriftSmartWriter<ferrum::thrift::TJSONProtocol>* tsw,
 const Vector1D& usage
 );
