#ifndef FERRUM_LIBNAR_SVI_UTIL_H_
#define FERRUM_LIBNAR_SVI_UTIL_H_

#include <vector>
#include <Eigen/Dense>

namespace ferrum {
  class StepSizeUpdater {
  private:
    double forgetting_;
    double delay_;
    double cached_;
    unsigned int iteration_;
    void _update();
  public:
    StepSizeUpdater();
    StepSizeUpdater(double f, double d);
    StepSizeUpdater& operator++();
    double operator()();
    double operator()() const;
    unsigned int iteration();
    unsigned int iteration() const;
    void forgetting(double f);
    void delay(double d);
  };

  void update_global_params
  (
   const Eigen::Matrix<double, 1, Eigen::Dynamic>& hypers,
   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat,
   int row,
   const Eigen::Matrix<double, 1, Eigen::Dynamic>& accumulated,
   double inter1,
   double inter2
   );
  void update_global_params
  (
   const std::vector<double>& hypers,
   std::vector<double>* old,
   std::vector<double>* accumulated,
   double inter1,
   double inter2
   );
  void update_global_params
  (
   std::vector<double>* old,
   std::vector<double>* accumulated,
   double inter1,
   double inter2
   );
}

#endif
