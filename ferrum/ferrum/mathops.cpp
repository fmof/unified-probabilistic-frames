#include "ferrum/mathops.hpp"

#ifdef USE_GSL_EXP
#include <gsl/gsl_sf_exp.h>
#else
#include <cmath>
#endif

namespace mathops {
  const double NEGATIVE_INFINITY = - std::numeric_limits<double>::infinity();
  const gsl_rng_type *which_gsl_rng = gsl_rng_mt19937;
  gsl_rng *rnd_gen = gsl_rng_alloc(which_gsl_rng);

  const double exp(const double x) {
    double res = 
#ifdef USE_GSL_EXP
      gsl_sf_exp(x);
#else
    std::exp(x);
#endif
    return res;
  }
  const double log(const double x) {
    double res = std::log(x);
    return res;
  }
  double log_sum_exp(const double* log_weights, int size) {
    double max = -std::numeric_limits<double>::infinity();
    for(int i = 0; i < size; ++i) {
      if(log_weights[i] > max)
	max = log_weights[i];
    }
    std::vector<double> tmp(size);
    double sum = 0.0;
    for(int i = 0; i < size; ++i) {
      double t = log_weights[i] - max;
      t = mathops::exp(t);
      sum += t;
    }
    return max + mathops::log(sum);
  }
}

namespace ferrum {

  ThreadRng::ThreadRng() {
    r_ = gsl_rng_alloc(gsl_rng_mt19937);
  }
  ThreadRng::~ThreadRng() {
    gsl_rng_free(r_);
  }
  gsl_rng* ThreadRng::get() {
    return r_;
  }

  namespace functor {
    double log(double x) {
      return mathops::log(x);
    }
    double exp(double x) {
      return mathops::exp(x);
    }

    Gaussian1D::Gaussian1D(const gsl_rng* rand_gen) :
      rg_(rand_gen),
      sigma_(1.0) {
    }
  }

  void prob_from_unnorm_lp
  (
   std::vector<double>* ulp
   ) {
    double lnorm = mathops::log_sum_exp(*ulp);
    ferrum::sum(-1 * lnorm, ulp);
    ferrum::exp(ulp);
  }
}
