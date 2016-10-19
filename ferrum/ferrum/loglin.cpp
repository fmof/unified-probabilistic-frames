#include "ferrum/loglin.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/mathops.hpp"
#include <vector>

namespace loglin {
  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>::UnigramMaxent(bool renormalize_when_set) :
    stale_log_normalizer_(true),
    renormalize_when_set_(renormalize_when_set),
    weights_(NULL),
    log_normalizer_(0.0) {
  }
  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>::UnigramMaxent() :
    UnigramMaxent<SupportType, WeightType>(true) {
  }

  // copy constructor
  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>::UnigramMaxent
  (
   const UnigramMaxent<SupportType,WeightType>& other
   ) :
    stale_log_normalizer_(other.stale_log_normalizer_),
    renormalize_when_set_(other.renormalize_when_set_),
    weights_(other.weights_),
    log_normalizer_(other.log_normalizer_) {
    //INFO << __PRETTY_FUNCTION__  << ": this.weights_ " << weights_ << ", other.weights_ " << other.weights_;
  }
  // move constructor
  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>::UnigramMaxent
  (
   UnigramMaxent<SupportType,WeightType>&& other
   ) :
    stale_log_normalizer_(other.stale_log_normalizer_),
    renormalize_when_set_(other.renormalize_when_set_),
    weights_(other.weights_),
    log_normalizer_(other.log_normalizer_) {
    //INFO << __PRETTY_FUNCTION__  << ": this.weights_ " << weights_ << ", other.weights_ " << other.weights_;
    other.weights_ = NULL;
  }

  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>&
  UnigramMaxent<SupportType, WeightType>::operator=
  (
   const UnigramMaxent<SupportType,WeightType>& other
   ){
    if(&other == this) return *this;
    stale_log_normalizer_ = other.stale_log_normalizer_;
    renormalize_when_set_ = other.renormalize_when_set_;
    log_normalizer_ = other.log_normalizer_;
    weights_ = other.weights_;
    return *this;
  }
  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>&
  UnigramMaxent<SupportType, WeightType>::operator=
  (
   UnigramMaxent<SupportType,WeightType>&& other
   ){
    if(&other == this) return *this;
    stale_log_normalizer_ = other.stale_log_normalizer_;
    renormalize_when_set_ = other.renormalize_when_set_;
    log_normalizer_ = other.log_normalizer_;
    weights_ = other.weights_;
    other.weights_ = NULL;
    return *this;
  }

  template <typename SupportType, typename WeightType>
  UnigramMaxent<SupportType, WeightType>::~UnigramMaxent() {
    // no-op
  }

  template <typename SupportType, typename WeightType>
  WeightType* UnigramMaxent<SupportType, WeightType>::weights() {
    return weights_;
  }
  template <typename SupportType, typename WeightType>
  void UnigramMaxent<SupportType, WeightType>::weights(WeightType* n_weights) {
    weights_ = n_weights;
    stale_log_normalizer_ = true;
    if(renormalize_when_set_) this->renormalize();
  }
  template <typename SupportType, typename WeightType>
  void UnigramMaxent<SupportType, WeightType>::weights(WeightType& n_weights) {
    weights_ = &n_weights;
    stale_log_normalizer_ = true;
    if(renormalize_when_set_) this->renormalize();
  }
  template <typename SupportType, typename WeightType>
  bool UnigramMaxent<SupportType, WeightType>::renormalize_when_set() {
    return renormalize_when_set_;
  }

  template <typename SupportType, typename WeightType>
  void UnigramMaxent<SupportType, WeightType>::renormalize() {
    log_normalizer_ = mathops::log_sum_exp(*weights_);
    stale_log_normalizer_ = false;
  }
  template <typename SupportType, typename WeightType>
  void UnigramMaxent<SupportType, WeightType>::renormalize_with_Z(double Z) {
    log_normalizer_ = Z;
    stale_log_normalizer_ = false;
  }

  template <typename SupportType, typename WeightType>
  double UnigramMaxent<SupportType, WeightType>::log_normalizer() {
    if(stale_log_normalizer_) {
      renormalize();
    }
    return log_normalizer_;
  }

  template <typename SupportType, typename WeightType> 
  double UnigramMaxent<SupportType, WeightType>::lp(SupportType obj) {
    return weights_->operator[](obj) - log_normalizer();
  }
  template <typename SupportType, typename WeightType> 
  double UnigramMaxent<SupportType, WeightType>::p(SupportType obj) {
    return mathops::exp(lp(obj));
  }


  template <typename SupportType, typename WeightType> 
  template <typename SparseDataSet>
  double UnigramMaxent<SupportType, WeightType>::ll(const SparseDataSet& data) {
    renormalize();
    double ll = 0.0;
    for(const auto& obj : data) {
      ll += obj.second * lp(obj.first);
    }
    return ll;
  }

  template <typename SupportType, typename WeightType> 
  template <typename SparseDataSet>
  WeightType UnigramMaxent<SupportType, WeightType>::ll_grad(const SparseDataSet& data) {
    renormalize();
    WeightType grad;
    // note that we reference grad.end() in order to always have amortized O(1) insertion
    typename WeightType::iterator g_it = grad.end();
    double sum = 0.0;
    for(const auto& obj : data) {
      sum += obj.second;
    }
    for(const auto& obj : data) {
      double val = obj.second - (sum * p(obj.first));
      // these two separate statements are needed because initially, g_it == 0x0 (due to empty container)
      g_it = grad.emplace(g_it, val);
      // and update to point to the end (for O(1) insertion)
      ++g_it;
    }
    return grad;
  }

  template <typename SupportType, typename WeightType> 
  double UnigramMaxent<SupportType, WeightType>::ll_dense_data(const std::vector<double>& data, bool renorm) {
    if(renorm) {
      renormalize();
    }
    double ll = 0.0;
    const size_t data_size = data.size();
    for(size_t i = 0; i < data_size; ++i) {
      //DEBUG << "obj " << i << ", data[i] = " << data[i] << ", lp(i) = " << lp(i);
      ll += data[i] * lp(i);
    }
    return ll;
  }
  template <typename SupportType, typename WeightType> 
  double UnigramMaxent<SupportType, WeightType>::ll_sparse_data(const std::vector<int>& counts, const std::vector<SupportType>& items, bool renorm) {
    if(renorm) {
      renormalize();
    }
    double ll = 0.0;
    //assert(counts.size() == items.size());
    const size_t data_size = items.size();
    for(size_t i = 0; i < data_size; ++i) {
      //DEBUG << "obj " << i << ", data[i] = " << data[i] << ", lp(i) = " << lp(i);
      ll += counts[i] * lp( items[i] );
    }
    return ll;
  }

  template <typename SupportType, typename WeightType> 
  WeightType UnigramMaxent<SupportType, WeightType>::ll_grad_dense_data(const std::vector<double>& data) {
    renormalize();
    WeightType grad;
    // note that we reference grad.end() in order to always have amortized O(1) insertion
    typename WeightType::iterator g_it = grad.end();
    const size_t data_size = data.size();
    const double N = ferrum::sum(data);
    DEBUG << "Total number of instances = " << N;
    for(size_t i = 0; i < data_size; ++i) {
      double val = data[i] - N*p(i);
      // these two separate statements are needed because initially, g_it == 0x0 (due to empty container)
      g_it = grad.emplace(g_it, val);
      // and update to point to the end (for O(1) insertion)
      ++g_it;
    }
    return grad;
  }

  template <typename SupportType, typename WeightType> 
  void UnigramMaxent<SupportType, WeightType>::ll_grad_dense_data
  (
   const std::vector<double>& data,
   double* grad   
   ) {
    renormalize();
    const size_t data_size = data.size();
    // note that we reference grad.end() in order to always have amortized O(1) insertion
    //typename WeightType::iterator g_it = grad.end();
    const double N = ferrum::sum(data);
    DEBUG << "Total number of instances = " << N;
    for(size_t i = 0; i < data_size; ++i) {
      double val = data[i] - N*p(i);
      grad[i] = val;
      // // these two separate statements are needed because initially, g_it == 0x0 (due to empty container)
      // g_it = grad.emplace(g_it, val);
      // // and update to point to the end (for O(1) insertion)
      // ++g_it;
    }
  }

  template class UnigramMaxent<int, std::vector<double> >;
  template class UnigramMaxent<size_t, std::vector<double> >;

  template double UnigramMaxent<size_t, std::vector<double> >::ll(const std::map<size_t, int>&);
  template double UnigramMaxent<int, std::vector<double> >::ll(const std::map<int, int>&);
  template std::vector<double> UnigramMaxent<size_t, std::vector<double> >::ll_grad(const std::map<size_t, int>&);
  template std::vector<double> UnigramMaxent<int, std::vector<double> >::ll_grad(const std::map<int, int>&);
}
