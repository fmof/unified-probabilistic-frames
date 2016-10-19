#include "ferrum/mathops.hpp"
#include "ferrum/stats.hpp"

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <Eigen/Dense>

namespace ferrum {
  double gamma::entropy(double a, double b) {
    double r = gsl_sf_psi(a);
    r *= (1.0 - a);
    r += a;
    r += gsl_sf_log(b);
    r += gsl_sf_lngamma(a);
    return r;
  }
  double gamma::entropy(std::pair<double, double> p) {
    return gamma::entropy(p.first, p.second);
  }

  ////////////////////////////////////////////////////////
  //// static function
  ////////////////////////////////////////////////////////

  ExpFamSampler::ExpFamSampler(bool mt) :
    mt_(mt) {
  }
  ExpFamSampler::ExpFamSampler() :
    ExpFamSampler(false) {
  }

  std::vector<double>
  ExpFamSampler::multinomial_uniform_dist(double uweight, int support_size) {
    return std::vector<double>(support_size, uweight/(double)support_size);
  }
  std::vector<double>
  ExpFamSampler::discrete_pc(double weight, const std::vector<double>& hyper) {
    return ferrum::sum(hyper,
		       std::vector<double>(hyper.size(),
					   weight/(double)hyper.size()));
  }
  std::vector<double>
  ExpFamSampler::multinomial_mean_uniform_pc(const std::vector<double>& hyper, double N) {
    if(mt_) {
      ERROR << __func__ << " is being called in a multithreaded context, but that isn't implemented. You may get weird results.";
    }
    const double d = (double)N/(double)hyper.size();
    std::vector<double> vec(hyper.size(), d);
    ferrum::sum_in_first(&vec, hyper);
    mathops::add_uniform_noise(&vec, -d, d);
    return vec;
  }

  std::vector<double>
  ExpFamSampler::multinomial_mean_uniform_dist(const std::vector<double>& hyper, double N) {
    if(mt_) {
      ERROR << __func__ << " is being called in a multithreaded context, but that isn't implemented. You may get weird results.";
    }
    const double d = (double)N/(double)hyper.size();
    std::vector<double> vec(hyper.size(), d);
    ferrum::sum_in_first(&vec, hyper);
    mathops::add_uniform_noise(&vec, -d, d);
    double norm = ferrum::sum(vec);
    ferrum::scalar_product(1.0/norm, &vec);
    return vec;
  }
  std::vector<double>
  ExpFamSampler::multinomial_mean_uniform_dist(size_t dim, double N) {
    if(mt_) {
      ERROR << __func__ << " is being called in a multithreaded context, but that isn't implemented. You may get weird results.";
    }
    const double d = (double)N/(double)dim;
    std::vector<double> vec(dim, d);
    mathops::add_uniform_noise(&vec, -d, d);
    double norm = ferrum::sum(vec);
    ferrum::scalar_product(1.0/norm, &vec);
    return vec;
  }

  // sampling code from http://lost-found-wandering.blogspot.com/2011/05/sampling-from-multivariate-normal-in-c.html and http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
  void
  MultiNormal::sample(int num_samples,
		      const Eigen::VectorXd& mean,
		      const MultiNormal::EMat& covar,
		      const gsl_rng* rng,
		      MultiNormal::EMat& samples) {
    ferrum::functor::Gaussian1D gauss_generator(rng);
    Eigen::LLT<EMat> cholSolver(covar);
    EMat transform;
    if (cholSolver.info()==Eigen::Success) {
      transform = cholSolver.matrixL();
    } else {
      WARN << "Failed computing the Cholesky decomposition. Falling back to eigen solver";
      Eigen::SelfAdjointEigenSolver<EMat> eigen_solver(covar);
      transform =
	eigen_solver.eigenvectors()*
	eigen_solver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
    }
    samples =
      (transform *
       EMat::NullaryExpr(mean.rows(),
			 num_samples,
			 gauss_generator)).colwise() +
      mean;
  }
}
