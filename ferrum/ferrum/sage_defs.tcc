#ifndef PUBLIC_UPF
#ifndef CRTLDA_SAGE_DEFS_TCC_
#define CRTLDA_SAGE_DEFS_TCC_

#include "ferrum/cblas_cpp.hpp"
#include "ferrum/loglin.hpp"
#include "ferrum/optimize.hpp"
#include "ferrum/sage_defs.hpp"
#include "ferrum/svi_util.hpp"

#include <mutex>
#include <thread>

namespace ferrum {
  template <typename EtaType>
  template<class Archive>
  void SageTopic<EtaType>::serialize(Archive& ar, const unsigned int version) {
    ar & support_size_;
    ar & background_;
    ar & log_partition_;
    ar & eta_;
    ar & regularization_type_;
  } 
  
  template <typename EtaType>
  SageTopic<EtaType>::SageTopic() : 
		       background_(),
		       eta_(new EtaType),
		       summed_weights_(new std::vector<double>()),
		       maxent_(),
		       tau_hyper_(1.0),
		       log_partition_(0.0),
		       support_size_(0),
		       regularization_type_(SageTopicRegularization::L2),
		       opt_how_(optimize::OptimizationMethod::LBFGS),
		       //minimizer_(std::shared_ptr<optimize::LibLBFGSMinimizer>(new optimize::LibLBFGSMinimizer)),
		       minimizer_(std::shared_ptr<optimize::LibLBFGSMinimizer>(NULL)),
		       sparsity_thres_(0.0),
		       sparse_(false) {
  }
  
  template <typename EtaType>
  SageTopic<EtaType>::SageTopic
  (
   int support_size,
   double kappa, 
   SageTopicRegularization reg_type
   ) :
    SageTopic<EtaType>(support_size,
		       kappa,
		       reg_type,
		       optimize::OptimizationMethod::LBFGS) {
  }
  template <typename EtaType>
  SageTopic<EtaType>::SageTopic
  (
   int support_size,
   double kappa, 
   SageTopicRegularization reg_type,
   optimize::OptimizationMethod opt_method
   ) : 
    background_(),
    eta_(new EtaType(support_size, 0.0)),
    //tau_(TauType(support_size, 1.0)),
    summed_weights_(new std::vector<double>(support_size, 0.0)),
    maxent_(),
    tau_hyper_(kappa),
    log_partition_(0.0),
    support_size_(support_size), 
    regularization_type_(reg_type),
    opt_how_(opt_method),
    sparsity_thres_(0.0),
    sparse_(false) {
    maxent_.weights(eta_.get());
    minimizer_ = optimize::Minimizer::make(opt_how_, support_size_);
  }
  template <typename EtaType>
  SageTopic<EtaType>::~SageTopic() {
    //INFO << __PRETTY_FUNCTION__ << " Deleting " << this << " :: " << summed_weights_;
    delete summed_weights_;
  }
  template <typename EtaType>
  SageTopic<EtaType>::SageTopic(const SageTopic<EtaType>& other) :
    background_( other.background_),
    eta_( new Eta( *(other.eta_) ) ),
    //tau_(other.tau_),
    summed_weights_( new std::vector<double>( *(other.summed_weights_) ) ),
    maxent_( other.maxent_ ),
    tau_hyper_(other.tau_hyper_),
    log_partition_( other.log_partition_ ),
    support_size_( other.support_size_),
    regularization_type_( other.regularization_type_ ),
    opt_how_( other.opt_how_ ),
    minimizer_( other.minimizer_ ),
    sparsity_thres_(0.0),
    sparse_(false) {
    //INFO << __PRETTY_FUNCTION__ << this << " from " << &other << " :: "  << ": this.summed_weights_ " << summed_weights_ << ", other.summed_weights_ " << other.summed_weights_;
    maxent_.weights(summed_weights_);
  }
  template <typename EtaType>
  SageTopic<EtaType>::SageTopic(SageTopic<EtaType>&& other) : 
    background_( std::move(other.background_) ),
    eta_( std::move(other.eta_) ),
    //tau_( std::move(other.tau_) ),
    summed_weights_ ( other.summed_weights_ ),
    maxent_ ( std::move( other.maxent_ ) ),
    tau_hyper_(other.tau_hyper_),
    log_partition_(other.log_partition_),
    support_size_(other.support_size_),
    regularization_type_(std::move(other.regularization_type_)),
    opt_how_(std::move(other.opt_how_)),
    minimizer_(std::move(other.minimizer_)),
    sparsity_thres_(other.sparsity_thres_),
    sparse_(other.sparse_) {
    //INFO << __PRETTY_FUNCTION__ << this << " from " << &other << " :: " << ": this.summed_weights_ " << summed_weights_ << ", other.summed_weights_ " << other.summed_weights_;
    other.summed_weights_ = NULL;
  }
  template <typename EtaType>
  SageTopic<EtaType>&
  SageTopic<EtaType>::operator=(const SageTopic<EtaType>& other) {
    background_ = other.background_;
    *eta_ = *(other.eta_);
    //tau_ = other.tau_;
    *summed_weights_ = *(other.summed_weights_);
    maxent_ = other.maxent_ ;
    maxent_.weights( summed_weights_ );
    tau_hyper_ = other.tau_hyper_;
    log_partition_ = other.log_partition_;
    support_size_ =  other.support_size_;
    regularization_type_ = other.regularization_type_;
    opt_how_ = other.opt_how_;
    minimizer_ = other.minimizer_;
    sparsity_thres_ = other.sparsity_thres_;
    sparse_ = other.sparse_;
    return *this;
  }
  template <typename EtaType>
  SageTopic<EtaType>&
  SageTopic<EtaType>::operator=(SageTopic<EtaType>&& other) {
    if(this != &other) {
      background_ = std::move(other.background_);
      eta_ = std::move(other.eta_);
      //tau_ = std::move(other.tau_);
      summed_weights_ = other.summed_weights_;
      maxent_ = std::move(other.maxent_);
      tau_hyper_ = other.tau_hyper_;
      log_partition_ = other.log_partition_;
      support_size_ =  other.support_size_;
      regularization_type_ = std::move(other.regularization_type_);
      opt_how_ = std::move(other.opt_how_);
      minimizer_ = std::move(other.minimizer_);
      sparsity_thres_ = other.sparsity_thres_;
      sparse_ = other.sparse_;
      other.summed_weights_ = NULL;
    }
    return *this;
  }
  template <typename EtaType>
  double& SageTopic<EtaType>::operator[](size_t idx) {
    return eta_->operator[](idx);
  }
  template <typename EtaType>
  void SageTopic<EtaType>::prepare() {
    if(summed_weights_->size() == 0) {
      summed_weights_->resize(support_size_, 0.0);
      if(minimizer_ == NULL) {
	minimizer_ = optimize::Minimizer::make(opt_how_, support_size_);
      }
    }
  }
  template <typename EtaType>
  void SageTopic<EtaType>::sparse(bool b) {
    sparse_ = b;
  }
  template <typename EtaType>
  void SageTopic<EtaType>::sparsity_threshold(double d) {
    sparsity_thres_ = d;
  }
  template <typename EtaType>
  bool SageTopic<EtaType>::sparse() {
    return sparse_;
  }
  template <typename EtaType>
  double SageTopic<EtaType>::sparsity_threshold() {
    return sparsity_thres_;
  }
  template <typename EtaType>
  minsky::Frame SageTopic<EtaType>::create_minsky_frame(int which_voc) {
    minsky::Frame f;
    minsky::Distribution d;
    d.__set_support_size(support_size_);
    minsky::Weights w;
    if(sparse_) {
      for(size_t i = 0; i < (size_t)support_size_; ++i) {
	double val = eta_->operator[](i);
	if(std::abs(val) >= sparsity_thres_) {
	  w.residual.push_back(val);
	  d.items.push_back((int)i);
	}
	d.__isset.items = true;
      }
    } else {
      std::copy(eta_->begin(), eta_->end(),
		std::back_inserter(w.residual));
    }
    w.__isset.residual = true;
    d.__set_sparse(sparse_);
    d.__set_sparsity_threshold(sparsity_thres_);
    d.__set_weights(w);
    d.__set_vocab_idx(which_voc);
    f.__set_distr(d);
    return f;
  }
  template <typename EtaType>
  void SageTopic<EtaType>::push_to_maxent() {
    this->eta(*eta_);
  }
  template <typename EtaType>
  typename SageTopic<EtaType>::MaxentModel*
  SageTopic<EtaType>::maxent_model() {
    return &maxent_;
  }
  template <typename EtaType>
  std::vector<double>* SageTopic<EtaType>::background() {
    return background_.get();
  }
  template <typename EtaType>
  void SageTopic<EtaType>::background(const std::shared_ptr<std::vector<double> >& background) {
    background_ = background;
  }
  template <typename EtaType>
  double SageTopic<EtaType>::log_partition() {
    return log_partition_;
  }
  // template <typename EtaType>
  // void SageTopic<EtaType>::renormalize(const std::vector<double>& logterm_sums) {
  //   log_partition_ = mathops::log_sum_exp(logterm_sums);
  // }
  // template <typename EtaType>
  // void SageTopic<EtaType>::_renormalize(bool skip_sum) {
  //   //summed_weights_ = ferrum::sum(*background_, *eta_);
  //   if(! skip_sum) {
  //     summed_weights_->assign(background_->begin(), background_->end());
  //     ferrum::sum_in_first(summed_weights_.get(), *eta_);
  //   }
  //   renormalize(*summed_weights_);
  // }
  template <typename EtaType>
  void SageTopic<EtaType>::renormalize() {
    summed_weights_->assign(background_->begin(), background_->end());
    ferrum::sum_in_first(summed_weights_, *eta_);
    log_partition_ = mathops::log_sum_exp(*summed_weights_);
    if(std::abs(log_partition_ - maxent_.log_normalizer()) > 1E-10) {
      maxent_.weights(summed_weights_);
      if(! maxent_.renormalize_when_set() ) {
	maxent_.renormalize_with_Z(log_partition_);
      }
    }
  }
  template <typename EtaType>
  void SageTopic<EtaType>::__eta() {
    const size_t size = eta_->size();
    for(size_t i = 0; i < size; ++i) {
      summed_weights_->operator[](i) = 
	background_->operator[](i) + eta_->operator[](i);
    }
    maxent_.weights(summed_weights_);
    log_partition_ = mathops::log_sum_exp(*summed_weights_);
    if(! maxent_.renormalize_when_set() ) {
      maxent_.renormalize_with_Z(log_partition_);
    }
  }
  template <typename EtaType>
  void SageTopic<EtaType>::eta(const EtaType& eta) {
    eta_->assign(eta.begin(), eta.end());
    this->__eta();
  }
  template <typename EtaType>
  EtaType SageTopic<EtaType>::eta() {
    return *eta_;
  }
  template <typename EtaType>
  double SageTopic<EtaType>::l_probability(const MaxentSupportType& obj) {
    return maxent_.lp(obj);
  }
  template <typename EtaType>
  double SageTopic<EtaType>::probability(const MaxentSupportType& obj) {
    return maxent_.p(obj);
  }

  template <typename EtaType>
  template <typename OutputType>
  OutputType SageTopic<EtaType>::as() {
    std::vector<double> logterm_sums = ferrum::sum(*background_, *eta_);
    renormalize(logterm_sums);
    ferrum::sum(-1.0 * log_partition_, &logterm_sums);
    ferrum::exp(&logterm_sums);
    OutputType out(logterm_sums);
    return out;
  }

  template <typename EtaType>
  template <typename OutputType>
  OutputType SageTopic<EtaType>::eta_as(bool normalize) const {
    if(normalize) {
      double lp = mathops::log_sum_exp(*eta_);
      OutputType res = ferrum::sum(-1.0 * lp, *eta_);
      ferrum::exp(&res);
      return res;
    } else {
      OutputType res(*eta_);
      return res;
    }
  }
  template <typename EtaType>
  template <typename OutputType>
  OutputType SageTopic<EtaType>::eta_as(bool normalize) {
    return const_cast< const SageTopic<EtaType>* >(this)->eta_as< OutputType >(normalize);
  }

  template <typename EtaType>
  double SageTopic<EtaType>::elbo_static(ClosureType* closure,
					 const double* weights,
					 const int n,
					 bool renorm) {
    double ll = 0.0;
    if(closure->topic_counts != NULL) {
      ll = closure->maxent->ll_dense_data(*(closure->topic_counts), renorm);
    } else if(closure->sparse_which_items != NULL && closure->sparse_which_counts != NULL) {
      ll = closure->maxent->ll_sparse_data(*(closure->sparse_which_counts), *(closure->sparse_which_counts), renorm);
    } else {
      ERROR << "Could not operate on closure";
      throw 5;
    }
    if(closure->prev_eta == NULL) {
      ERROR << "Could not operate on closure with null prev_eta";
      throw 6;
    }
    // now compute the regularizer
    const double reg_multiplier = closure->regularizer_multiplier;
    double weight_sum = 0.0;
    for(int i = 0; i < n; ++i) {
      const auto val = weights[i];
      auto rval = compute_sage_regularizer(val, closure->regularizer_type, closure->prev_eta->operator[](i));
      weight_sum += rval;
    }
    double regularizer_s = reg_multiplier * weight_sum;
    ll += (-.5 * regularizer_s);
    return -ll;
  }

  template <typename EtaType>
  double SageTopic<EtaType>::elbo(const std::vector<double>& counts, bool renorm) {
    ClosureType ct;
    ct.topic_counts = const_cast<std::vector<double>*>(&counts);
    ct.background = this->background();
    ct.regularizer_multiplier = 1.0;
    ct.maxent = this->maxent_model();
    ct.regularizer_type = regularization_type_;
    double ll = SageTopic<EtaType>::elbo_static(&ct, eta_.get()->data(), eta_->size(), renorm);
    return ll;
  }
  template <typename EtaType>
  double SageTopic<EtaType>::elbo(const std::vector<int>& counts, const std::vector<int>& which, bool renorm) {
    ClosureType ct;
    ct.topic_counts = NULL;
    ct.sparse_which_counts = const_cast<std::vector<int>*>(&counts);
    ct.sparse_which_items = const_cast<std::vector<int>*>(&which);
    ct.background = this->background();
    ct.regularizer_multiplier = 1.0;
    ct.maxent = this->maxent_model();
    ct.regularizer_type = regularization_type_;
    double ll = SageTopic<EtaType>::elbo_static(&ct, eta_.get()->data(), eta_->size(), renorm);
    return ll;
  }


  /**
   * static function
   */
  template <typename EtaType>
  double SageTopic<EtaType>::elbo_eval(void *fparams,
				       const double* point,
				       const int size) {
    ClosureType* closure = (ClosureType*)fparams;
    // update the maxent model
    // add the trial_weights to the background
    std::vector<double> nweights = 
      optimize::LibLBFGSVector::sum(point, closure->background, size);
    // reset the maxent weights; this will renormalize everything
    closure->maxent->weights(nweights);
    double ll = elbo_static(closure, nweights.data(), size,
			    !(closure->maxent->renormalize_when_set()));
    return ll;
  }

  /**
   * static function
   */
  template <typename EtaType>
  void SageTopic<EtaType>::elbo_grad(void *fparams,
				     const double* point,
				     double* grad,
				     const int size) {
    ClosureType* closure = (ClosureType*)fparams;
    std::vector<double> nweights;
    if(! (closure->is_model_ready) ) {
      nweights.resize(size, 0.0);
      cblas_dcopy(size, point, 1, nweights.data(), 1);
      cblas_daxpy(size, 1.0, closure->background->data(), 1,
		  nweights.data(), 1);
      // reset the maxent weights; this will renormalize everything
      closure->maxent->weights(nweights);
      closure->is_model_ready = true;
    }
    catlas_dset(size, 0.0, grad, 1);
    closure->maxent->ll_grad_dense_data(*(closure->topic_counts), grad);
    // subtract off regularization
    const double reg_multiplier = closure->regularizer_multiplier * 0.5;
    for(int i = 0; i < size; ++i) {
      const double val = point[i];
      const auto gval = compute_grad_sage_regularizer(val, closure->regularizer_type,
						      closure->prev_eta->operator[](i));
      grad[i] = -grad[i];
      grad[i] += (reg_multiplier * gval);
      //lbfgs_grad[i] = -grad[i] + reg_multiplier * gval;
    }
  }

  // this is a static function
  template <typename EtaType>
  double SageTopic<EtaType>::liblbfgs_elbo(void *fparams, const lbfgsfloatval_t *trial_weights,
					   lbfgsfloatval_t *lbfgs_grad, const int n, const lbfgsfloatval_t step) {
    ClosureType* closure = (ClosureType*)fparams;
    // update the maxent model
    // add the trial_weights to the background
    std::vector<double> nweights = 
      optimize::LibLBFGSVector::sum(trial_weights, closure->background, n);
    // reset the maxent weights; this will renormalize everything
    closure->maxent->weights(nweights);
    double ll = elbo_static(closure, nweights.data(), n, !(closure->maxent->renormalize_when_set()));

    closure->is_model_ready = true;
    elbo_grad(fparams, trial_weights, lbfgs_grad, n);
    closure->is_model_ready = false;
    return ll;
  }
  
  template <typename EtaType>
  std::shared_ptr<optimize::Function>
  SageTopic<EtaType>::get_optimizable_func(ClosureType* params) {
    std::shared_ptr<optimize::Function> func;
    if(opt_how_ == optimize::OptimizationMethod::LBFGS) {
      func = std::shared_ptr<optimize::Function>
	(new optimize::LibLBFGSFunction(&SageTopic<EtaType>::liblbfgs_elbo,
					&optimize::LibLBFGSNoOp::progress,
					(void*)params));
    }
    else {
      switch(opt_how_) {
      case optimize::OptimizationMethod::ADAGRAD:
	{
	  optimize::MinimizationFunctionCreator<optimize::OptimizationMethod::ADAGRAD> mfc;
	  func = mfc(&SageTopic<EtaType>::elbo_eval,
		     &SageTopic<EtaType>::elbo_grad,
		     (void*)params);
	}
	break;
      default:
	WARN << "Trying to create function for optimization method " << opt_how_ << "; unexpected, but it might work";
	break;
      }
      // func = optimize::Minimizer::get_function(opt_how_,
      // 					       &SageTopic<EtaType>::elbo_grad,
      // 					       (void*)params);
    }
    return func;
  }

  // TODO: the following two methods could be refactored a bit
  template <typename EtaType>
  int SageTopic<EtaType>::fit_topic(const std::vector<double>* counts, double reg_mult) {
    ClosureType params;
    params.topic_counts = const_cast< std::vector<double>* >(counts);
    params.background = this->background();
    params.regularizer_multiplier = reg_mult;
    params.maxent = this->maxent_model();
    params.regularizer_type = regularization_type_;
    Eta prev_eta(this->eta());
    params.prev_eta = &prev_eta;
    //optimize::LibLBFGSFunction my_func = this->get_liblbfgs_func(&params);
    std::shared_ptr<optimize::Function> my_func = get_optimizable_func(&params);
    Eta point = this->get_optimization_initial_point();
    //optimize::LibLBFGSMinimizer optimizer(point.size());
    //int opt_status = optimizer.minimize(&my_func, point);
    int opt_status = minimizer_->minimize(my_func.get(), point);
    // set the inner eta to the point (i.e., do not call te private method)
    this->eta(point);
    return opt_status;
  }
  template <typename EtaType>
  int SageTopic<EtaType>::fit_topic_svi(const std::vector<double>* counts, double reg_mult, double i1, double i2) {
    if(opt_how_ == optimize::OptimizationMethod::LBFGS) {
      ClosureType params;
      params.topic_counts = const_cast< std::vector<double>* >(counts);
      params.background = this->background();
      params.regularizer_multiplier = reg_mult;
      params.maxent = this->maxent_model();
      params.regularizer_type = regularization_type_;
      Eta prev_eta(this->eta());
      params.prev_eta = &prev_eta;
      //optimize::LibLBFGSFunction my_func = this->get_liblbfgs_func(&params);
      std::shared_ptr<optimize::Function> func = get_optimizable_func(&params);
      Eta point = this->get_optimization_initial_point();
      //optimize::LibLBFGSMinimizer optimizer(point.size());
      //int opt_status = optimizer.minimize(&my_func, point);
      int opt_status = minimizer_->minimize(func.get(), point);
      // eta_ = i1 * eta_ + i2 * point
      update_global_params
	(
	 eta_.get(), // the old point
	 &point, // once this function returns, point is now all 0
	 i1,
	 i2
	 );
      // call the private method, since we've already set the inner
      // eta_ values
      this->__eta();
      return opt_status;
    } else if(opt_how_ == optimize::OptimizationMethod::ADAGRAD) {
    } else {
      ERROR << "Unknown optimization method " << opt_how_;
      throw 10;
    }
    return 0;
  }

  template <typename EtaType>
  double SageTopic<EtaType>::variance_first_deriv(const std::pair<double, double>& pair, double e_i) {
    double deriv = (.5 + pair.first) * gsl_sf_psi(pair.first);
    double asq = (pair.first - 1) * (pair.first - 1);
    if(asq < 1E-8) {
      asq = 1E-8;
      //WARN << "Capping (a-1)**2 to 1E-8";
    }
    double esq = e_i * e_i;
    if(esq < 1E-10) {
      esq = 1E-10;
      //WARN < "Capping eta**2 to 1E-10";
    }
    deriv += (.5 * esq / pair.second / asq);
    deriv -= 1;
    return deriv;
  }
  template <typename EtaType>
  double SageTopic<EtaType>::variance_second_deriv(const std::pair<double, double>& pair, double e_i) {
    double deriv = (.5 + pair.first) * gsl_sf_psi_n(2, pair.first);
    double am1 = pair.first - 1;
    double acu = am1 * am1 * am1;
    if(acu < 1E-8) {
      acu = 1E-8;
      //WARN << "Capping (a-1)**3 to 1E-8";
    }
    double esq = e_i * e_i;
    if(esq < 1E-10) {
      esq = 1E-10;
      //WARN < "Capping eta**2 to 1E-10";
    }
    deriv += (esq / pair.second / acu);
    return deriv;
  }

  template <typename EtaType>
  void SageTopic<EtaType>::update_variances(std::vector<std::pair<double,double> >* params, double i1, double i2) {
    size_t i = 0;
    for(std::pair<double,double>& pair : *params) {
      const double e_i = eta_->operator[](i);
      double n_a = pair.first;
      double alpha = 0.2; // as recommended from original SAGE code
      for(int iter = 0; iter < 10; ++iter) {
	double deriv1 = variance_first_deriv(pair, e_i);
	double deriv2 = variance_second_deriv(pair, e_i);
	double deriv = deriv1/deriv2;
	if(std::isnan(deriv)) {
	  WARN << "Optimization of variance unstable, got NaN; ending optimization";
	  break;
	}
	if(std::isinf(deriv)) {
	  WARN << "Optimization of variance unstable, got " << deriv << "; ending optimization";
	  break;
	}
	pair.first -= (alpha * deriv1/deriv2);
      }
      if(i1 == 0.0) {
	pair.first = n_a;
	pair.second = e_i * e_i / (pair.first - 1);
      } else {
	// pair.first = i1*old_value + i2*newly_computed_value
	pair.first *= i2;
	pair.first += (i1 * n_a);
	pair.second *= i2;
	pair.second += (i1 * (e_i * e_i / (pair.first - 1)) );
      }
      ++i;
    }
  }

  // This gets an initial point.
  // Following the original SAGE implementation, this initializes to the zero vector
  // (which makes sense, because in expectation, due to the sparsity-inducing prior, 
  // it should be zero).
  template <typename EtaType>
  typename SageTopic<EtaType>::Eta
  SageTopic<EtaType>::get_optimization_initial_point() {
    return Eta(support_size_, 0.0);
  }
  template <typename EtaType>
  int SageTopic<EtaType>::support_size() {
    return support_size_;
  }

  template <typename CorpusT>
  std::vector<double> SageInitializer::topic_fit_subset
  (
   const std::vector<double>& hyper,
   CorpusT* corpus,
   std::shared_ptr< std::vector<double> > background,
   minsky::MinskyAnnotationWrapper maw
   // minsky::AnnotationLevel::type annot_level,
   // minsky::StructureType::type st_level
   ) {
    const double avg_words_per_doc = (double)num_words_ / (double)(corpus->num_docs());
    int init_num = num_docs_for_init_;
    if(init_num <= 0) {
      init_num = (int)(10000.0 / avg_words_per_doc);
    }
    if(init_num > corpus->num_docs()) {
      WARN << "The number of desired documents for subset initialization (" << init_num << ") is >= the total number of documents available (" << corpus->num_docs() << "); setting to " << corpus->num_docs();
      init_num = corpus->num_docs();
    }
    // pick init_num docs
    std::vector<double> dense_counts(hyper.size(), 0.0);
    CorpusT* rand_corp = NULL;
    {
      std::lock_guard<std::mutex> lock(*num_calls.subset_mut);
      rand_corp = corpus->random_subset((size_t)init_num, ++(num_calls.topic_subset));
    }
    INFO << "Performing initialization on random subset " << rand_corp->get_name() << " of " << (rand_corp->num_docs()) << " documents";
    for(const auto& di : *rand_corp) {
      const auto& doc = *(di.document);
      for(const auto& pair : ferrum::doc_multinomial( doc, maw) ) {
	dense_counts[pair.first] += pair.second;
      }
    }
    delete rand_corp;
    typedef SageTopic<std::vector<double> > TopicType;
    TopicType topic(hyper.size(), 1.0, reg_topic_how_);
    topic.background(background);
    int opt_status = topic.fit_topic(&dense_counts, 1.0);
    INFO << "Maxent optimization in initialization of topic resulted in " << opt_status << " status";
    std::vector<double> point = topic.eta_as<std::vector<double> >(false);
    return point;
  }

  /**
   * Provide parameters eta over V elements for a SageTopic maxent model. 
   */
  template <typename CorpusT, typename Vocab>
  std::vector<double> SageInitializer::topic
  (
   const std::vector<double>& hyper,
   CorpusT* corpus,
   const Vocab* vocab,
   std::shared_ptr<std::vector<double> > background,
   minsky::MinskyAnnotationWrapper maw
   // minsky::AnnotationLevel::type annot_level,
   // minsky::StructureType::type st_level
   ) {
    switch(fit_topic_how_) {
    case TopicInitializerChoice::UNIFORM:
      return topic_uniform(hyper, background, (double)(vocab->num_words()));
    case TopicInitializerChoice::SUBSET:
      return topic_fit_subset(hyper, corpus, background, maw);
    case TopicInitializerChoice::BACKGROUND:
      return topic_background(hyper, background, (double)(vocab->num_words()));
    default:
      ERROR << "Invalid topic fit selection \"" << fit_topic_how_ << "\"";
      throw 4;
    }
  }  
}

#endif
#endif
