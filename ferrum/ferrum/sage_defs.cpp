#ifndef PUBLIC_UPF
#include "ferrum/sage_defs.hpp"

namespace ferrum {
  std::istream& operator>>(std::istream& in, SageTopicRegularization &how) {
    std::string token;
    in >> token;
    if (token == "L2")
      how = SageTopicRegularization::L2;
    else if (token == "IMPROPER" )
      how = SageTopicRegularization::IMPROPER;
    else {
      ERROR << "Unknown regularization choice " << token;
      throw 5;
    }
    //else throw boost::program_options::validation_error("Invalid unit");
    return in;
  }

  double compute_sage_regularizer(double value, SageTopicRegularization how, double prev_val) {
    double res = 0.0;
    switch(how) {
    case SageTopicRegularization::L2:
      res = value * value;
      break;
    case SageTopicRegularization::IMPROPER:
      res = value * prev_val * prev_val * value;
      break;
    default:
      ERROR << "Unknown regularization type " << how;
      throw 4;
    }
    return res;
  }
  double compute_grad_sage_regularizer(double value, SageTopicRegularization how, double prev_val) {
    double res = 0.0;
    switch(how) {
    case SageTopicRegularization::L2:
      res = 2.0 * value;
      break;
    case SageTopicRegularization::IMPROPER:
      res = 2.0 * value * prev_val;
      break;
    default:
      ERROR << "Unknown regularization type " << how;
      throw 4;
    }
    return res;
  }

  std::istream& operator>>(std::istream& in, TopicInitializerChoice &how) {
    std::string token;
    in >> token;
    if (token == "UNIFORM")
      how = TopicInitializerChoice::UNIFORM;
    else if (token == "SUBSET")
      how = TopicInitializerChoice::SUBSET;
    else if (token == "BACKGROUND")
      how = TopicInitializerChoice::BACKGROUND;
    else {
      ERROR << "Unknown initialization choice " << token;
      throw 5;
    }
    //else throw boost::program_options::validation_error("Invalid unit");
    return in;
  }
  /////////////////////////////////////////////////////////////

  SageInitializer::SageInitializer
  (
   int num_words,
   int num_threads,
   TopicInitializerChoice fit_how, 
   SageTopicRegularization reg_topic_how
   ) : 
    num_words_(num_words),
    num_threads_(num_threads),
    fit_topic_how_(fit_how),
    reg_topic_how_(reg_topic_how),
    num_docs_for_init_(-1)  {
    num_calls.topic_subset = 0;
    num_calls.subset_mut = std::shared_ptr<std::mutex>(new std::mutex);
  }
  SageInitializer::SageInitializer() :
    SageInitializer(-1, // num words
		    1, // num threads
		    TopicInitializerChoice::SUBSET, 
		    SageTopicRegularization::L2) {
  }
  SageInitializer::SageInitializer(int num_words) :
    SageInitializer(num_words, 1,
		    TopicInitializerChoice::SUBSET,
		    SageTopicRegularization::L2) {
  }
  SageInitializer::SageInitializer(int num_words, int num_threads) :
    SageInitializer(num_words, num_threads,
		    TopicInitializerChoice::SUBSET,
		    SageTopicRegularization::L2) {
  }
  int SageInitializer::num_threads(bool safe) {
    if(safe) {
      return (num_threads_ <= 0) ? 1 : num_threads_;
    }
    return num_threads_;
  }
  void SageInitializer::num_docs_for_init(int nd) {
    num_docs_for_init_ = nd;
  }
  void SageInitializer::sparse(bool thrift, bool opt, double thresh) {
    sparsity_info.serialize_ = thrift;
    sparsity_info.opt_ = opt;
    sparsity_info.threshold_ = thresh;
  }
  bool SageInitializer::sparse() {
    return sparsity_info.opt_;
  }
  double SageInitializer::threshold() {
    return sparsity_info.threshold_;
  }
  std::vector<std::pair<double, double> > SageInitializer::gamma_variance(size_t size) {
    return std::vector<std::pair<double, double> >(size, std::make_pair<double, double>(1.0, 1.0));
  }
  std::vector<std::pair<double, double> > SageInitializer::verb_variance(size_t size) {
    return gamma_variance(size);
  }
  std::vector<std::pair<double, double> > SageInitializer::arc_variance(size_t size) {
    return gamma_variance(size);
  }
  const SageTopicRegularization& SageInitializer::regularize_how() const {
    return reg_topic_how_;
  }

  std::vector<double> SageInitializer::usage(const std::vector<double>& hyper, const double N) {
    const double d = (double)N/(double)hyper.size();
    std::vector<double> vec(hyper.size(), d);
    ferrum::sum_in_first(&vec, hyper);
    mathops::add_uniform_noise(&vec, -d, d);
    return vec;
  }
  std::vector<double> SageInitializer::usage_templates(const std::vector<double>& hyper, const double N) {
    return usage(hyper, N);
  }
  std::vector<double> SageInitializer::usage(const std::vector<double>& hyper) {
    if(num_words_ <= 0) {
      ERROR << "Attempting to call .usage with num_words in corpus == " << num_words_;
      throw 3;
    }
    return usage(hyper, (double)num_words_);
  }
  std::vector<double> SageInitializer::slots(const std::vector<double>& hyper) {
    if(num_words_ <= 0) {
      ERROR << "Attempting to call " << __func__ << " with num_words in corpus == " << num_words_;
      throw 3;
    }
    return usage(hyper, (double)num_words_);
  }

  std::vector<double> SageInitializer::assignment(size_t num) {
    double inv_nt_ = 1.0/(double)num;
    std::vector<double> vec(num, inv_nt_);
    mathops::add_uniform_noise(&vec, -inv_nt_, inv_nt_);
    double norm = ferrum::sum(vec);
    ferrum::scalar_product(1.0/norm, &vec);
    return vec;
  }


  std::vector<double> SageInitializer::topic_uniform
  (
   const std::vector<double>& hyper,
   std::shared_ptr< std::vector<double> > background,
   double num_types
   ) {
    const double iv = 1.0/num_types;
    std::vector<double> vec( *background );
    ferrum::scalar_product(-1.0, &vec);
    mathops::add_uniform_noise(&vec, -iv, iv);
    return vec;
  }

  std::vector<double> SageInitializer::topic_background
  (
   const std::vector<double>& hyper,
   std::shared_ptr<std::vector<double> > background,
   double num_types
   ) {
    const double iv = 1.0/num_types;
    std::vector<double> vec( (size_t)num_types, 0.0 );
    mathops::add_uniform_noise(&vec, -iv, iv);
    return vec;
  }

  void complete_background_init
  (
   int num_words,
   std::vector<double>* raw_bck,
   double min_log_prob,
   double min_freq
   ) {
    // make sure there are no zero counts
    ferrum::ensure_min(min_freq, raw_bck);
    // divide by the total number of words in the corpus
    ferrum::scalar_product(1.0/(double)num_words, raw_bck);
    // now log it all
    ferrum::log(raw_bck);
    // truncate to min_log_prob
    ferrum::ensure_min(min_log_prob, raw_bck);
    // renormalize
    dmc::cat::log_renormalize(raw_bck);
  }

}

template class ferrum::SageTopic< std::vector<double> >;
#endif
