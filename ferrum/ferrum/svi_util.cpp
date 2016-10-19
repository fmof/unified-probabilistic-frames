#include "ferrum/svi_util.hpp"
#include "ferrum/util.hpp"

#include <cmath>
#include <Eigen/Dense>
#include <vector>

namespace ferrum {
  void update_global_params
  (
   const Eigen::Matrix<double, 1, Eigen::Dynamic>& hypers,
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat,
   int row_index,
   const Eigen::Matrix<double, 1, Eigen::Dynamic>& accumulated,
   double inter1,
   double inter2
   ) {
    if(inter1 == 0.0) {
      mat.row(row_index) = hypers + accumulated;
    } else {
      mat.row(row_index) += inter1 * mat.row(row_index) +
	inter2 * (hypers + accumulated);
    }
  }
  void update_global_params
  (
   const std::vector<double>& hypers,
   std::vector<double>* old,
   std::vector<double>* accumulated,
   double inter1,
   double inter2
   ) {
    if(inter1 == 0.0) {
      *old = hypers;
      ferrum::sum_in_first(old, *accumulated);
    } else {
      ferrum::sum_in_first(accumulated, hypers);
      ferrum::linear_combination_in_first
	(
	 old,
	 *accumulated,
	 inter1,
	 inter2
	 );
    }
    ferrum::set(accumulated, 0.0);
  }
  void update_global_params
  (
   std::vector<double>* old,
   std::vector<double>* accumulated,
   double inter1,
   double inter2
   ) {
    if(inter1 == 0.0) {
      std::copy(accumulated->begin(), accumulated->end(),
		old->begin());
    } else {
      ferrum::set(old, 0.0);
      ferrum::linear_combination_in_first
	(
	 old,
	 *accumulated,
	 inter1,
	 inter2
	 );
    }
    ferrum::set(accumulated, 0.0);
  }
  StepSizeUpdater::StepSizeUpdater(double f, double d) :
    forgetting_(f),
    delay_(d),
    cached_(0.0),
    iteration_(0.0) {
    _update();
  }
  StepSizeUpdater::StepSizeUpdater() : StepSizeUpdater(1.0, 1.0) {
  }
  StepSizeUpdater& StepSizeUpdater::operator++() {
    ++iteration_;
    _update();
    return *this;
  }
  void StepSizeUpdater::_update() {
    const double base = (iteration_ + delay_);
    const double eb = std::pow(base, -forgetting_);
    cached_ = eb;
  }
  double StepSizeUpdater::operator()() {
    return cached_;
  }
  double StepSizeUpdater::operator()() const {
    return cached_;
  }
  unsigned int StepSizeUpdater::iteration() {
    return iteration_;
  }
  unsigned int StepSizeUpdater::iteration() const {
    return iteration_;
  }
  void StepSizeUpdater::forgetting(double f) {
    forgetting_ = f;
  }
  void StepSizeUpdater::delay(double d) {
    delay_ = d;
  }
}
