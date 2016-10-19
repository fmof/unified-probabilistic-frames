#ifndef FERRUM_LOGLIN_H_ 
#define FERRUM_LOGLIN_H_

#include "ferrum/logging.hpp"
#include "ferrum/optimize.hpp"
#include "ferrum/util.hpp"
#include <vector>

namespace loglin {
  /**
   * Note that WeightType must be indexable by SupportType
   */
  template <typename SupportType, typename WeightType>
  class UnigramMaxent {
  private:
    bool stale_log_normalizer_;
    bool renormalize_when_set_;
  protected:
    WeightType* weights_; // this really should be a std::shared_ptr
    double log_normalizer_;
  public:
    UnigramMaxent();
    UnigramMaxent(bool renormalize_when_set);
    ~UnigramMaxent();
    UnigramMaxent(const UnigramMaxent<SupportType,WeightType>& other);
    UnigramMaxent(UnigramMaxent<SupportType,WeightType>&& other);
    UnigramMaxent& operator=(const UnigramMaxent<SupportType,WeightType>& other);
    UnigramMaxent& operator=(UnigramMaxent<SupportType,WeightType>&& other);
    double lp(SupportType obj);
    double p(SupportType obj);
    double log_normalizer();

    template <typename SparseDataSet> double ll(const SparseDataSet& data);
    template <typename SparseDataSet> WeightType ll_grad(const SparseDataSet& data);
    double ll_dense_data(const std::vector<double>& data, bool renormalize = true);
    double ll_sparse_data(const std::vector<int>& data, const std::vector<SupportType>& items, bool renormalize = true);
    WeightType ll_grad_dense_data(const std::vector<double>& data);
    void ll_grad_dense_data(const std::vector<double>& data, double* grad);
    void weights(WeightType& n_weights);
    void weights(WeightType* n_weights);
    bool renormalize_when_set();
    WeightType* weights();
    void renormalize();
    void renormalize_with_Z(double Z);
  };

  struct IntUnigramMaxentClosure {
    typedef loglin::UnigramMaxent<int, std::vector<double> > MaxentModel;
    MaxentModel* model;
    std::vector<double>* counts;
    double regularizer_strength = 0.0;
  };
}

#endif
