#ifndef FERRUM_STATS_HPP_
#define FERRUM_STATS_HPP_

#include <utility>

#include <Eigen/Dense>
#include <gsl/gsl_rng.h>

namespace ferrum {
  /**
   * A Gamma distribution parametrized with shape and scale (not shape and rate)
   */
  class gamma {
  public:
    static double entropy(double a, double b);
    static double entropy(std::pair<double,double>);
  };

  class ExpFamSampler {
  private:
    bool mt_;
  public:
    ExpFamSampler();
    ExpFamSampler(bool mt);
    std::vector<double> multinomial_uniform_dist(double uweight, int support_size);
    std::vector<double> discrete_pc(double weight, const std::vector<double>& hyper);
    /**
     * Given hyper of size K, return
     * hyper + vec(N/K) + {Uniform(-N/K, N/K)}_k
     */
    std::vector<double> multinomial_mean_uniform_pc(const std::vector<double>& hyper, double N);
    /**
     * Given hyper of size K, return a *normalized* form of
     * hyper + vec(N/K) + {Uniform(-N/K, N/K)}_k
     */
    std::vector<double> multinomial_mean_uniform_dist(const std::vector<double>& hyper, double N);
    /**
     * Given hyper of size K, return a *normalized* form of
     * vec(N/K) + {Uniform(-N/K, N/K)}_k
     */
    std::vector<double> multinomial_mean_uniform_dist(size_t dim, double N);
  };

  class MultiNormal {
  public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMat;

    /**
     * Sample num_samples from the specified multivariate Gaussian.
     * mean is a *column* vector of (K x 1), covar is the K x K covariance matrix.
     * samples will be reset to a K x num_samples matrix
     */
    static void sample(int num_samples, const Eigen::VectorXd& mean, const EMat& covar, const gsl_rng* rg, EMat& samples);
  };
}

#endif
