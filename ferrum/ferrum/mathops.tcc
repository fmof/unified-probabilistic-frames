#ifndef FERRUM_MATHOPS_TCC_
#define FERRUM_MATHOPS_TCC_

#include "ferrum/util.hpp"

#include <vector>

namespace mathops {
  template <typename Container>
  inline const double log_add(const Container& log_probs) {
    double sum = mathops::NEGATIVE_INFINITY;
    for(double lp : log_probs) {
      sum = mathops::log_add(sum, lp);
    }
    return sum;
  }
  inline const double log_add(const std::vector<double>& log_probs) {
    double sum = mathops::NEGATIVE_INFINITY;
    const size_t si = log_probs.size();
    const double* const lp_ptr = log_probs.data();
    for(size_t i = 0; i < si; ++i) {
      sum = mathops::log_add(sum, *(lp_ptr + i));
    }
    return sum;
  }

  template <typename Container>
  inline const double log_sum_exp(const Container& log_weights) {
    typedef typename Container::value_type T;
    const T& max = ferrum::max(log_weights);
    Container nterms = ferrum::sum(-1 * max, log_weights);
    return (double)max + log_add(nterms);
  }

  template <>
  inline const double log_sum_exp(const std::vector<double>& log_weights) {
    const double max = ferrum::max(log_weights);
    const double* const dp = log_weights.data();
    const size_t size = log_weights.size();
    std::vector<double> tmp(size);
    double* const tp = tmp.data();
    for(size_t i = 0; i < size; ++i) {
      double t = dp[i] - max;
      tp[i] = exp(t);
    }
    double sum = ferrum::sum(tmp);
    return max + log(sum);
  }
}

namespace ferrum {
  namespace functor {
    template <typename Index>
    const double Gaussian1D::operator()(Index, Index) const {
      return gsl_ran_gaussian(rg_, sigma_);
    }
  }
}

#endif
