#ifndef FERRUM_LIBNAR_CRTLDA_DMC_VARIATIONAL_H_
#define FERRUM_LIBNAR_CRTLDA_DMC_VARIATIONAL_H_

#include "ferrum/concrete.hpp"
#include "concrete_util/uuid_util.h"
#include "ferrum/crtlda_variational.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/dmc.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/svi_util.hpp"
#include "ferrum/util.hpp"

//#include "ferrum/minsky.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <map>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ferrum {
  /**
   * Variational inference for the discrete kinds template model.
   */
  template <typename DType>
  class DiscreteCRTLDAVariational {
  private:
    typedef UniformHyperSeedWeightedInitializer DiscreteVariationalInitializer;
    typedef std::vector<double> TopicType;
    typedef std::string GV;
    typedef std::string RV;
    typedef std::vector<double> ModelTopicType;
    typedef DiscreteModelGlobalSlots M;
    typedef std::vector< double > Vector1D;
    typedef std::vector< std::vector< double > > Vector2D;
    int num_docs_;

    Vector2D var_template_usage_params_;
    // GLOBAL: Templates x Slots
    Vector2D var_template_slot_params_;
    // GLOBAL: Templates x Verbs
    Vector2D var_template_verb_params_;
    // GLOBAL: Slots x arcs
    Vector2D var_slot_arc_params_;

    // GLOBAL:
    //Vector1D buffer_template_usage_params_;
    // GLOBAL: Templates x Slots
    Vector2D buffer_template_slot_params_;
    // GLOBAL: Templates x Verbs
    Vector2D buffer_template_verb_params_;
    // GLOBAL: Slots x Arcs
    Vector2D buffer_slot_arc_params_;

    // the model itself
    M* model_;

    int num_templates_;
    int num_slots_;
    //int num_frames_;
    //int num_roles_;
    int num_verbs_;
    int num_arcs_;

    Vector1D template_usage_hypers_;
    Vector1D slot_usage_hypers_;
    //Vector1D frame_hypers_;
    //Vector1D role_hypers_;
    Vector1D verb_hypers_;
    Vector1D arc_hypers_;

    // These are sparse, repeated views. There may be some value x for which, e.g.,
    // frame_type_counts[0][0][i] == frame_type_counts[0][0][j] == x. 
    // All ${x}s_in_docs_.size() == ${x}_type_counts_.size(), for every dimension.
    // Note that if x occurs X times in a given entity e, then
    // sum_{i : frames_in_docs_[0][e][i] == x} frame_type_counts[0][e][i] == X.
    // Moreover, note that because of the dependence of verbs on frames (arcs on roles),
    // the frame/verb (role/arc) lists are kept in sync. That is, 
    // frames_in_docs_[][].size() == verbs_in_docs_[][].size()
    // However, the frame/role (verb/arc) lists are not necessarily in sync:
    // frames_in_docs_[][].size() != roles_in_docs_[][].size()
    // Of course, sum_i frame_type_counts_[][][i] == sum_i role_type_counts_[][][i]
    std::vector< size_t > num_ents_in_docs_;
    std::vector< std::vector< std::vector< int > > > verbs_in_docs_;
    std::vector< std::vector< std::vector< int > > > verb_type_counts_;
    std::vector< std::vector< std::vector< int > > > arcs_in_docs_;
    std::vector< std::vector< std::vector< int > > > arc_type_counts_;

    DiscreteVariationalInitializer vi_init_;

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ERROR << "This isn't complete";
      throw 5;
    } 

  public:
    DiscreteCRTLDAVariational<DType>() :  num_docs_(0), model_(NULL) {
    }
    DiscreteCRTLDAVariational<DType>(M* model) :
    num_docs_(0),
	       model_(model),
	       num_templates_(model->num_templates()), num_slots_(model->num_slots()),
      template_usage_hypers_(model->hyper_theta()), 
      slot_usage_hypers_(model->hyper_slot()),
		verb_hypers_(model->hyper_gov()), arc_hypers_(model->hyper_rel()) {
    }
    ~DiscreteCRTLDAVariational() {
    }

    void num_templates(int nt) {
      num_templates_ = nt;
    }
    int num_templates() {
      return num_templates_;
    }
    void num_slots(int ns) {
      num_slots_ = ns;
    }
    int num_slots() {
      return num_slots_;
    }

    void model(M* model) {
      model_ = model;
    }
    M reconstruct_model() {
      M model(this->num_topics_);
      ERROR << "This method is not complete";
      throw 3;
      // model.hyper_theta(this->usage_hypers_);
      // model.hyper_word(this->word_hypers_);
      // this->update_model(&model);
      return model;
    }

    template <typename WCountType = int, typename RecordType = WCountType>
    RecordType populate_obs_lists(const std::map<int, WCountType>& doc_multi,
    				  std::vector< std::vector<int> >* words_in_docs,
    				  std::vector< std::vector<RecordType> >* word_type_counts) {
      std::vector<WCountType> word_types;
      std::vector<int> wcounts;
      RecordType num_words = (RecordType)0.0;
      for(const auto& count_pair : doc_multi) {
    	const int word = count_pair.first;
    	word_types.push_back(word);
    	RecordType c = (RecordType)(count_pair.second);
    	wcounts.push_back(c);
    	num_words += c;
      }
      words_in_docs->push_back(word_types);
      word_type_counts->push_back(wcounts);
      return num_words;
    }

    /**
     * DiscreteVariationalInitializer must have four functions:
     */
    void init(const DiscreteVariationalInitializer& vi,
	      const Vocabulary<GV>& verb_vocab,
	      const Vocabulary<RV>& arc_vocab) {
      this->vi_init_ = vi;
      num_verbs_ = verb_vocab.num_words();
      num_arcs_ = arc_vocab.num_words();
      //buffer_template_usage_params_.resize(num_templates_, 0.0);
      for(int t = 0; t < num_templates_; ++t) {
	var_template_slot_params_.push_back(vi_init_.slots(model_->hyper_slot()));
	var_template_verb_params_.push_back(vi_init_.verbs(model_->hyper_verb()));
	buffer_template_slot_params_.push_back(Vector1D(num_slots_));
	buffer_template_verb_params_.push_back(Vector1D(num_verbs_));
      }
      for(int s = 0; s < num_slots_; ++s) {
	var_slot_arc_params_.push_back(vi_init_.arcs(model_->hyper_arc()));
	buffer_slot_arc_params_.push_back(Vector1D(num_arcs_));
      }
    }

    // inline void get_usage_estimates(std::vector<std::vector<double> >*  usage_ptr) {
    //   for(const auto& use : var_template_usage_params_) {
    // 	const double norm = ferrum::sum(use);
    // 	usage_ptr->push_back(ferrum::scalar_product(1.0/norm, use));
    //   }
    // }

    Vector2D get_distributions(const Vector2D* ptr) {
      Vector2D topics;
      typedef Vector2D::value_type Inner;
      for(const auto& utopic : *ptr) {
	const double norm = ferrum::sum(utopic);
	Inner topic(ferrum::scalar_product(1.0/norm, utopic));
	topics.push_back(topic);
      }
      return topics;
    }

    void update_model(M* model, bool heldout = false) {
      if(!heldout) {
	model->prior_gov(get_distributions(&var_template_verb_params_));
	model->prior_slot_usage(get_distributions(&var_template_slot_params_));
	model->prior_rel(get_distributions(&var_slot_arc_params_));
      }
      //model->prior_template_usage(get_distributions(&var_template_usage_params_));
    }
    void update_model(bool heldout = false) {
      update_model(model_, heldout);
    }
    // transfer the learned parameters back to the model
    void transfer_learned_parameters(M* model) {
      update_model(model, true);
    }
    void transfer_learned_parameters() {
      transfer_learned_parameters(model_);
    }
    void update_template_assignments
    (
     Vector1D* thematic_params,
     const Vector1D& narrative_params,
     const Vector1D& expected_thematic_counts,
     const Vector1D& expected_syntactic_counts
     ) {
      ferrum::copy(&narrative_params, thematic_params);
      for(size_t idx = 0; idx < num_templates_; ++idx) {
       	thematic_params->operator[](idx) += 
       	  expected_thematic_counts[idx] +
	  expected_syntactic_counts[idx];
      }
      ferrum::prob_from_unnorm_lp(thematic_params);
    }
    void update_slot_assignments
    (
     Vector1D* thematic_params,
     const Vector1D& expected_thematic_counts,
     const Vector1D& expected_syntactic_counts
     ) {
      ferrum::copy(&expected_thematic_counts, thematic_params);
      ferrum::sum_in_first(thematic_params, expected_syntactic_counts);
      ferrum::prob_from_unnorm_lp(thematic_params);
    }

    void update_template_usage(std::vector<double>* var_assign,
			       const std::vector<double>& buffer) {
      *var_assign = template_usage_hypers_;
      ferrum::sum_in_first(var_assign, buffer);
    }

    void _update_template_slot_usage(int index_template, double i1 = 0.0, double i2 = 1.0) {
      update_global_params
	(
	 model_->hyper_slot(),
	 &(var_template_slot_params_[index_template]),
	 &(buffer_template_slot_params_[index_template]),
	 i1,
	 i2
	 );
    }
    void _update_template_verb_usage(int index_template, double i1 = 0.0, double i2 = 1.0) {
      update_global_params
	(
	 model_->hyper_verb(),
	 &(var_template_verb_params_[index_template]),
	 &(buffer_template_verb_params_[index_template]),
	 i1,
	 i2
	 );
    }
    void _update_slot_arc_usage(int index_slot, double i1 = 0.0, double i2 = 1.0) {
      update_global_params
	(
	 model_->hyper_arc(),
	 &(var_slot_arc_params_[index_slot]),
	 &(buffer_slot_arc_params_[index_slot]),
	 i1,
	 i2
	 );
    }

    void update_template_slot_usage(int num_threads, double i1 = 0.0, double i2 = 1.0) {
      omp_set_num_threads(num_threads);
#pragma omp parallel for
      for(int ti = 0; ti < num_templates_; ++ti) {	  
	_update_template_slot_usage(ti);
      }
    }
    void update_template_verb_usage(int num_threads, double i1 = 0.0, double i2 = 1.0) {
      omp_set_num_threads(num_threads);
#pragma omp parallel for
      for(int ti = 0; ti < num_templates_; ++ti) {
	_update_template_verb_usage(ti);
      }
    }
    void update_slot_arc_usage(int num_threads, double i1 = 0.0, double i2 = 1.0) {
      omp_set_num_threads(num_threads);
#pragma omp parallel for
      for(int ti = 0; ti < num_slots_; ++ti) {	  
	_update_slot_arc_usage(ti);
      }
    }
    
    template <typename BSP>
    void write_situations
    (
     const DType& doc,
     BSP protocol_sp,
     const Vector2D& t_assign,
     const Vector2D& s_assign
     ) {
      // assert(doc.num_entities() == t_assign.size());
      // assert(t_assign.size() == s_assign.size());
      // // this is currently a "mode" serializer
      // Vector1D t_max = ferrum::row_arg_max(t_assign);
      // Vector1D s_max = ferrum::row_arg_max(s_assign);
      // const size_t ne = doc.num_entities();
      // concrete::util::uuid_factory uf;
      // for(size_t i = 0; i < ne; ++i) {
      // 	std::string situation_kind = std::to_string(t_max[i]);
      // 	std::string situation_role = std::to_string(s_max[i]);
      // 	// now iterate over the mentions
      // 	{
      // 	  //const auto& conc_mention;
      // 	  concrete::SituationMention sm;
      // 	  concrete::UUID uuid;
      // 	  uuid.__set_uuidString(uf.get_uuid());
      // 	  sm.__set_uuid(uuid);
      // 	  sm.__set_situationKind(situation_kind);
      // 	  std::vector<concrete::MentionArgument> sm_args;
      // 	  // fill
      // 	  concrete::MentionArgument ma;
      // 	  ma.__set_role(std::to_string(s_max[i]));
      // 	  concrete::UUID ma_ent_uuid;
      // 	  // make this better
      // 	  ma_ent_uuid.__set_uuidString(conc_mention.id);
      // 	  ma.__set_entityMentionId(ma_ent_uuid);
      // 	  sm_args.push_back(ma);
      // 	  sm.__set_argumentList(sm_args);
      // 	  // optional?
      // 	  sm.__set_situationType("EVENT");
      // 	}
      // }
    }

    template <typename SL> // SL == ferrum::SituationLabeler<>
    void e_step
    (
     const VStrategy& strategy,
     int learn_iter,
     const size_t batch_start,
     const size_t batch_end,
     SL& sit_lab
     ) {
      Vector2D slots_grads(num_templates_, Vector1D(num_slots_, 0.0));
      Vector2D verbs_grads(num_templates_, Vector1D(num_verbs_, 0.0));
      Vector2D arcs_grads(num_slots_, Vector1D(num_arcs_, 0.0));
      if(strategy.em_verbosity >= 1) {
	INFO << "E-step";
      }
      for(size_t ti = 0; ti < num_templates_; ++ti) {
	dmc::dirichlet::grad_log_partition_static(var_template_slot_params_[ti], &(slots_grads[ti]));
	dmc::dirichlet::grad_log_partition_static(var_template_verb_params_[ti], &(verbs_grads[ti]));
      }
      for(size_t si = 0; si < num_slots_; ++si) {
	dmc::dirichlet::grad_log_partition_static(var_slot_arc_params_[si], &(arcs_grads[si]));
      }
#ifdef _PRINT_CSV_
      std::cout << "FOO,e.iter,doc.idx,ent.idx,var,var.id,assign.val\n";
#endif
      const size_t batch_size = batch_end - batch_start;
      omp_set_num_threads(strategy.num_e_threads);
#pragma omp parallel for
      for(size_t di = 0; di < batch_size; ++di){
	const size_t num_ents = num_ents_in_docs_[di];
	Vector2D var_template_assign_params(num_ents);
	Vector2D var_slot_assign_params(num_ents);
	Vector2D expected_verb_counts(num_ents, Vector1D(num_templates_, 0.0));
	Vector2D expected_arc_counts(num_ents, Vector1D(num_slots_, 0.0));
	// grad(A(var{slot dists})) * s_{d,e}
	Vector1D expected_slot_for_templates(num_templates_, 0.0);
	// t_{d,e} * grad(A(var{slot dists})) ==
	// trans(grad(A(var{slot dists}))) * t_{d,e}
	Vector1D expected_template_for_slots(num_slots_, 0.0);
	Vector1D butp(num_templates_, 0.0);
	for(size_t ent_idx = 0; ent_idx < num_ents; ++ent_idx) {
	  var_template_assign_params[ent_idx] = vi_init_.assignment(num_templates_);
	  var_slot_assign_params[ent_idx] = vi_init_.assignment(num_slots_);
	  // update the innermost variables first
	  // use two separate loops to account for sparse view
	  const size_t num_ments = verbs_in_docs_[di][ent_idx].size();
	  for(size_t ment_idx = 0; ment_idx < num_ments; ++ment_idx) {
	    int verb_count = verb_type_counts_[di][ent_idx][ment_idx];
	    // accumulate expected counts, for the specific frame of interest (i.e., the one that was actually observed)
	    int obs_verb = verbs_in_docs_[di][ent_idx][ment_idx];
	    auto fg_col = ferrum::column(verbs_grads, obs_verb);
	    ferrum::linear_combination_in_first
	      (
	       &(expected_verb_counts[ent_idx]),
	       fg_col,
	       1.0,
	       (double)verb_count
	       );
	  }
	  for(size_t ment_idx = 0; ment_idx < num_ments; ++ment_idx) {
	    int arc_count = arc_type_counts_[di][ent_idx][ment_idx];
	    // accumulate expected counts, for the specific frame of interest (i.e., the one that was actually observed)
	    int obs_arc = arcs_in_docs_[di][ent_idx][ment_idx];
	    auto fg_col = ferrum::column(arcs_grads, obs_arc);
	    ferrum::linear_combination_in_first
	      (
	       &(expected_arc_counts[ent_idx]),
	       fg_col,
	       1.0,
	       (double)arc_count
	       );
	  }
	}
	for(int e_iter = 0; e_iter < strategy.num_e_iters; ++e_iter) {
	  if((strategy.em_verbosity >= 3) || 
	     (strategy.em_verbosity >= 2 && di % 1000 == 0)) {
	    INFO << "\tDocument " << di << " of iteration " << e_iter;
	  }
	  bool last_iter = ( (e_iter + 1) == strategy.num_e_iters);
	  if(strategy.em_verbosity >= 4) {
	    INFO << "\t\tE-step sub-iteration number " << e_iter;
	  }
	  std::vector<double> usage_grad = 
	    dmc::dirichlet::grad_log_partition_static(var_template_usage_params_[di]);
	  for(size_t ent_idx = 0; ent_idx < num_ents; ++ent_idx) {
#ifdef _PRINT_CSV_
	    if(e_iter == 0) {
	      for(size_t tii = 0; tii < num_templates_; ++tii) {
	    	std::cout << "FOO," << learn_iter << "." << -1 << "," << di << "," << ent_idx;
	    	std::cout << ",template," << tii << "," << var_template_assign_params[ent_idx][tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << -1 << "," << di << "," << ent_idx;
		std::cout << ",b.tusage," << tii << "," << buffer_template_usage_params_[tii] << "\n";
	    	std::cout << "FOO," << learn_iter << "." << -1 << "," << di << "," << ent_idx;
	    	std::cout << ",tusage," << tii << "," << var_template_usage_params_[di][tii] << "\n";
	    	std::cout << "FOO," << learn_iter << "." << -1 << "," << di << "," << ent_idx;
	    	std::cout << ",tusage.grad," << tii << "," << usage_grad[tii] << "\n";
	      }
	      for(size_t tii = 0; tii < num_slots_; ++tii) {
	    	std::cout << "FOO," << learn_iter << "." << -1 << "," << di << "," << ent_idx;
	    	std::cout << ",slot," << tii << "," << var_slot_assign_params[ent_idx][tii] << "\n";
	      }
	    }
#endif

	    ferrum::set(&expected_slot_for_templates, 0.0);
	    ferrum::product(slots_grads, var_slot_assign_params[ent_idx],
			     &expected_slot_for_templates);
	    update_template_assignments
	      (
	       &(var_template_assign_params[ent_idx]), 
	       usage_grad,
	       expected_slot_for_templates,
	       expected_verb_counts[ent_idx]
	       );
	    ferrum::sum_in_first(&butp, var_template_assign_params[ent_idx]);
	    // and accumulate into expected counts
	    ferrum::set(&expected_template_for_slots, 0.0);
	    ferrum::transpose_product(slots_grads,
				       var_template_assign_params[ent_idx],
				       &expected_template_for_slots);
	    update_slot_assignments
	      (
	       &(var_slot_assign_params[ent_idx]), 
	       expected_template_for_slots,
	       expected_arc_counts[ent_idx]
	       );
#ifdef _PRINT_CSV_
	    if(1) {
	      for(size_t tii = 0; tii < num_templates_; ++tii) {
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",template," << tii << "," << var_template_assign_params[ent_idx][tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",tusage.grad," << tii << "," << usage_grad[tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",b.tusage," << tii << "," << buffer_template_usage_params_[tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",tusage," << tii << "," << var_template_usage_params_[di][tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",esft," << tii << "," << expected_slot_for_templates[tii] << "\n";
	      }
	      for(size_t tii = 0; tii < num_slots_; ++tii) {
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",slot," << tii << "," << var_slot_assign_params[ent_idx][tii] << "\n";
		std::cout << "FOO," << learn_iter << "." << e_iter << "," << di << "," << ent_idx;
		std::cout << ",etfs," << tii << "," << expected_template_for_slots[tii] << "\n";
	      }
	    }
#endif

	    // update buffers if last time
	    if(last_iter) {
	      // update the buffers
	      {
		const size_t num_ments = verbs_in_docs_[di][ent_idx].size();
		for(size_t ti = 0; ti < num_templates_; ++ti) {
		  const double t_weight = var_template_assign_params[ent_idx][ti];
		  for(size_t mi = 0; mi < num_ments; ++mi) {
		    size_t obs_idx = verbs_in_docs_[di][ent_idx][mi];
		    int num_times = verb_type_counts_[di][ent_idx][mi];
		    buffer_template_verb_params_[ti][obs_idx] += num_times * t_weight;
		  }
		  for(size_t si = 0; si < num_slots_; ++si) {
		    ferrum::linear_combination_in_first(&(buffer_template_slot_params_[ti]),
							 var_slot_assign_params[ent_idx],
							 1.0,
							 t_weight);
		  }
		}
		for(size_t si = 0; si < num_slots_; ++si) {
		  const double s_weight = var_slot_assign_params[ent_idx][si];
		  for(size_t mi = 0; mi < num_ments; ++mi) {
		    size_t obs_idx = arcs_in_docs_[di][ent_idx][mi];
		    int num_times = arc_type_counts_[di][ent_idx][mi];
		    buffer_slot_arc_params_[si][obs_idx] += num_times * s_weight;
		  }
		}
	      }
	      if(sit_lab.do_labeling) {
		Corpus<DType>* corp_for_labeling =
		  static_cast<Corpus<DType>*>(sit_lab.unsafe_corpus);
		assert(corp_for_labeling != NULL);
		const auto& doc = corp_for_labeling->operator[](di);
		std::string suffix = "iter" + doc.id;
		// auto sit_writer = sit_lab.get(suffix); //boost::shared_ptr<P>
		// write_situations(doc, sit_writer, var_template_assign_params, var_slot_assign_params);
	      }
	    }
	  } // end for loop(ent_idx)
	  update_template_usage(&(var_template_usage_params_[di]), butp);
	  ferrum::set(&butp, 0.0);
	} // end for loop(e_iter)
      } // end for loop(doc \in corpus)
    } // end e_step()

    void m_step(VStrategy& strategy, unsigned int batch_size) {
      if(strategy.em_verbosity >= 1) {
	INFO << "M-step";
      }
      double rho = strategy.ssu();
      double divider = rho / (double)batch_size;
      update_template_slot_usage(strategy.num_m_threads, 1.0 - rho, divider);
      update_template_verb_usage(strategy.num_m_threads, 1.0 - rho, divider);
      update_slot_arc_usage(strategy.num_m_threads, 1.0 - rho, divider);
      ++(strategy.ssu);
    }

    // TODO: Fix this
    void update_hypers() {
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

    template <typename DKPrintingStruct>
    void print_dist_diag(const DKPrintingStruct* const print_struct,
			 const VStrategy& strategy,
			 int learn_iter, bool last_iter) {
      if(print_struct != NULL) {
	if(learn_iter % strategy.print_topics_every == 0 || last_iter) {
	  // if we haven't updated, then we should get some form of the current topics
	  auto* vocab_ptr = print_struct->g_vocab;
	  if(learn_iter % strategy.update_model_every) {
	    std::vector< ModelTopicType > topics;
	    for(const auto& utopic : var_template_verb_params_) {
	      const double norm = ferrum::sum(utopic);
	      TopicType topic(ferrum::scalar_product(1.0/norm, utopic));
	      topics.push_back(topic);
	    }
	    model_->print_verbs(strategy.print_topics_k, *vocab_ptr, topics);
	  } else {
	    model_->print_verbs(strategy.print_topics_k, *vocab_ptr);
	  }
	} else {
	  WARN << "The SmartWriter wrapper was not null, but the vocab wrapper was. We cannot print any vocab distributions (akin to topics) with out.";
	}
      } else {
	WARN << "Cannot print diagnostic template distributions to console because the vocabs are not provided.";
      }
    }

    template <typename DKPrintingStruct>
    void _print_in_learn
    (
     const VStrategy& strategy, 
     const DKPrintingStruct* const print_struct,
     DKVWriters* sw_wrapper,
     int epoch,
     int learn_iter,
     bool last_iter,
     bool model_changed
     ) {
      // TERMINAL PRINTING:: This is different than printing to file
      print_dist_diag(print_struct, strategy, learn_iter, last_iter);
      if(sw_wrapper == NULL) {
	return;
      }
      std::string suff = "epoch" + std::to_string(epoch) + ".emiter" + std::to_string(learn_iter);
      if(learn_iter % strategy.print_topics_every == 0 || last_iter) {
	// in order to print any of the distributions, we need the vocabs
	if(print_struct != NULL) {
	  // TODO: VOCAB DISTRIBUTION PRINTING
	  // TODO: make this bit more general (macro, perhaps?)
	  if(sw_wrapper->to_file_verb_usage()) {
	    ferrum::SmartWriter* swptr = sw_wrapper->sw_verb_usage();
	    std::ostream& out_stream = swptr->get(suff);
	    INFO << "Writing verb distributions to " << swptr->name();
	    // if we haven't updated, then we should get some form of the current topics
	    auto* vocab_ptr = print_struct->g_vocab;
	    if(learn_iter % strategy.update_model_every) {
	      std::vector< ModelTopicType > topics;
	      for(const auto& utopic : var_template_verb_params_) {
		const double norm = ferrum::sum(utopic);
		TopicType topic(ferrum::scalar_product(1.0/norm, utopic));
		topics.push_back(topic);
	      }
	      model_->print_verbs(num_verbs_, *vocab_ptr, out_stream, topics);
	    } else {
	      model_->print_verbs(num_verbs_, *vocab_ptr, out_stream);
	    }
	  }
	  if(sw_wrapper->to_model_tsv()) {
	    ferrum::SmartWriter* swptr = sw_wrapper->sw_model_tsv();
	    INFO << "Writing entire model as TSV to " << swptr->name();
	    std::ostream& out_stream = swptr->get(suff);
	    if(! (learn_iter % strategy.update_model_every)) {
	      update_model(strategy.heldout);
	    }
	    auto cp_ps = *print_struct;
	    cp_ps.num_per_gov = num_verbs_;
	    cp_ps.num_per_rel = num_arcs_;
	    model_->print_templates(out_stream, &cp_ps);
	  }
	} else {
	  WARN << "The SmartWriter wrapper was not null, but the vocab wrapper was. We cannot print any vocab distributions (akin to topics) with out.";
	}
      }
      // USAGE PRINTING
      if((learn_iter % strategy.print_usage_every == 0 && learn_iter > 0) || model_changed || last_iter) {
	if(sw_wrapper->to_file_template_usage()) {
	  ferrum::SmartWriter* swptr = sw_wrapper->sw_template_usage();
	  std::ostream& out_stream = swptr->get(suff);
	  INFO << "Writing document-template usage output to " << swptr->name();
	  //model_->print_template_usage(out_stream);
	  ferrum::print_2d_distribution(get_distributions(&var_template_usage_params_), out_stream);
	}
	if(sw_wrapper->to_file_slot_usage() ) {
	  ferrum::SmartWriter* swptr = sw_wrapper->sw_slot_usage();
	  std::ostream& out_stream = swptr->get(suff);
	  INFO << "Writing document-template usage output to " << swptr->name();
	  model_->print_slot_usage(out_stream);
	}
      }
    }

    template <typename DKPrintingStruct>
    void init_batch(const Corpus<DType>& corpus,
		    const DKPrintingStruct& print_struct)
    {
      num_docs_ = corpus.num_docs();
      // ensure that all document-specific globals are empty
      num_ents_in_docs_.resize(num_docs_, 0);
      var_template_usage_params_.clear();
      verbs_in_docs_.clear();
      verb_type_counts_.clear();
      arcs_in_docs_.clear();
      arc_type_counts_.clear();
      // and now initialize
      int num_words = 0;
      for(size_t di = 0; di < num_docs_; ++di) {
	const DType& doc = corpus[di];
	var_template_usage_params_.push_back(vi_init_.usage_template(model_->hyper_usage(), doc.num_entities()));
	num_ents_in_docs_[di] = doc.num_entities();
	for(const auto& entity : doc.entities()) {
	  verbs_in_docs_.push_back(std::vector<std::vector<int> >());
	  verb_type_counts_.push_back(std::vector<std::vector<int> >());
	  arcs_in_docs_.push_back(std::vector<std::vector<int> >());
	  arc_type_counts_.push_back(std::vector<std::vector<int> >());
	  num_words +=
	    populate_obs_lists(entity.gov_lemma_histogram(*print_struct.g_vocab),
			       &(verbs_in_docs_[di]),
			       &(verb_type_counts_[di]));
	  populate_obs_lists(entity.rel_str_histogram(*print_struct.r_vocab),
			     &(arcs_in_docs_[di]),
			     &(arc_type_counts_[di]));
	}
      }
    }

    void label
    (const size_t batch_start,
     const size_t batch_end
     ) {
      for(size_t di = batch_start; di < batch_end; ++di){
	const size_t num_ents = num_ents_in_docs_[di];
	for(size_t ent_idx = 0; ent_idx < num_ents; ++ent_idx) {
	}
      }
    }

    /**
     * DKPrintingStruct: e.g., ferrum::DiscreteKindPrinter
     */
    template <typename DKPrintingStruct>
    void learn(const Corpus<DType>& corpus,
	       VStrategy& strategy,
	       int epoch,
	       const DKPrintingStruct& print_struct,
	       DKVWriters* sw_wrapper = NULL) {
      bool last_iter = false;
      bool force_update = false;
      bool model_changed = false;
      // Initialize the batch
      init_batch(corpus, print_struct);
      if(strategy.em_verbosity >= 0) {
	INFO << "EM Epoch " << epoch;
      }
      const size_t batch_size = strategy.batch_size <= 0 ? num_docs_ : (size_t)strategy.batch_size;
      const size_t num_batches = 
	(size_t)((num_docs_ % batch_size) ?
		 (num_docs_ / batch_size + 1) :
		 (num_docs_ / batch_size));
      size_t corpus_start = 0, 
	corpus_end = batch_size;
      bool label_docs = (strategy.label_every > 0);
      SituationLabeler<concrete::util::TBinaryProtocol> sit_lab(label_docs);
      sit_lab.make_csw("situations.tbinary");
      if(label_docs) {
	sit_lab.unsafe_corpus =
	  static_cast<void*>
	  (
	   const_cast< Corpus<DType>* >(&corpus)
	   );
      }
      for(size_t batch_i = 0; batch_i < num_batches; ++batch_i) {
	e_step(strategy, epoch, corpus_start, corpus_end, sit_lab);
	if(! strategy.heldout) {
	  m_step(strategy, (unsigned int)batch_size);
	}	
	if( strategy.hyper_update_iter > 0 ) {//  && 
	  // ( (learn_iter > strategy.hyper_update_min 
	  //    && learn_iter % strategy.hyper_update_iter == 0 ) 
	  //   || last_iter) ) {
	  update_hypers();
	  force_update = model_changed = true;
	  //model_changed = true;
	}
	if(corpus_end % strategy.update_model_every == 0 || force_update || last_iter) {
	  update_model(strategy.heldout);
	  model_changed = true;
	}
	_print_in_learn(strategy, &print_struct, sw_wrapper, epoch, corpus_end, last_iter, model_changed);
	if(strategy.label_every > 0) {
	  //label();
	}
	update_model(strategy.heldout);
	corpus_start = corpus_end + 1;
	corpus_end += batch_size;
	if(corpus_end >= num_docs_) {
	  corpus_end = num_docs_;
	  last_iter = true;
	}
      }

    }

    void learn() {
      const VStrategy strategy;
      learn(strategy);
    }

  };
}

#endif
