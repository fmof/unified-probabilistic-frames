#include "ferrum/logging.hpp"
#include "ferrum/util.hpp"
#include "ferrum/version.hpp"

#include <fstream>
#include <iostream>
#include <ostream>
#include <memory>
#include <mutex>
#include <sys/types.h>
#include <unistd.h>

#include "ferrum/cblas_cpp.hpp"
#include <vector>

#include <algorithm>
#include <string>

namespace ferrum {

  std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec) {
    const size_t K = vec.size();
    out << "(";
    for(size_t k = 0; k < K; ++k) {
      out << vec.at(k);
      if(k + 1 < K) {
	out << " ";
      } else {
	out << ")";
      }
    }
    return out;
  }
  std::ostream& operator<<(std::ostream&& out, const std::vector<double>& vec) {
    const size_t K = vec.size();
    out << "(";
    for(size_t k = 0; k < K; ++k) {
      out << vec.at(k);
      if(k + 1 < K) {
	out << " ";
      } else {
	out << ")";
      }
    }
    return out;
  }

  void lower(std::string& data) {
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  }

  template <>
  void sum_in_first(std::vector< float >* const x, const std::vector< float >& y) {
    const size_t dim = x->size();
    if(dim != y.size()) {
      throw 1;
    }
    float* xp = x->data();
    cblas_saxpy(dim, 1, y.data(), 1, xp, 1);
  };

  template <>
  void sum_in_first(std::vector< double >* const x, const std::vector< double >& y) {
    const size_t dim = x->size();
    if(dim != y.size()) {
      throw 1;
    }
    double* xp = x->data();
    cblas_daxpy(dim, 1, y.data(), 1, xp, 1);
  };
  template <>
  void sum_in_first(std::vector< double >* const x, const Eigen::VectorXd& y) {
    const size_t dim = x->size();
    if(dim != (size_t)y.rows()) {
      throw 1;
    }
    double* xp = x->data();
    cblas_daxpy(dim, 1, y.data(), 1, xp, 1);
  };

  double variance(const Eigen::VectorXd& samples) {
    const double value_mean = samples.mean();
    const double num = (double)samples.size();
    double samp_var = (samples.array() - value_mean).square().sum()/num;
    return samp_var;
  }

  /**
   * Compute y = Ax
   */
  void product(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A,
	       const std::vector<double>& x,
	       std::vector<double>* y) {
    if(y->size() > 0 && y->size() != (size_t)A.rows()) {
      ERROR << "Cannot compute y = Ax, where y \\in R^" << (y->size()) << " and A \\in R^{" << A.rows() << ", " << A.cols() << "}";
      throw 5;
    }
    if(y->size() == 0) {
      WARN << "Passing in a result vector of length 0; resizing to " << A.rows();
      y->resize( A.rows() );
    }
    Eigen::VectorXd::Map(y->data(), y->size()) =
      A * Eigen::VectorXd::Map(x.data(), x.size());
  }

  /**
   * Compute y = Ax
   */
  void product(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A,
	       const double* x,
	       double* y) {
    Eigen::VectorXd::Map(y, A.rows()) =
      A * Eigen::VectorXd::Map(x, A.cols());
  }

  template <>
  void scalar_product(const double& x, std::vector< double >* y) {
    const size_t dim = y->size();
    cblas_dscal(dim, x, y->data(), 1);
  };
  template <>
  std::vector<double> scalar_product(const double& x, const std::vector< double >& y) {
    const size_t dim = y.size();
    std::vector<double> res(dim);
    cblas_dcopy(dim, y.data(), 1, res.data(), 1);
    cblas_dscal(dim, x, res.data(), 1);
    return res;
  };
  template <>
  void scalar_product(const int& x, std::vector< double >* y) {
    const size_t dim = y->size();
    cblas_dscal(dim, x, y->data(), 1);
  };
  template <>
  std::vector<double> scalar_product(const int& x, const std::vector< double >& y) {
    const size_t dim = y.size();
    std::vector<double> res(dim);
    cblas_dcopy(dim, y.data(), 1, res.data(), 1);
    cblas_dscal(dim, x, res.data(), 1);
    return res;
  };
  /**
   * ensure every value is >= fl
   * return the number of values changed (pre-normalization)
   */
  int floor_and_l1_norm_probs(std::vector<double>& vals, double fl) {
    // make sure there are no zero counts
    int num_changed = ferrum::ensure_min(fl, &vals);
    if(num_changed == 0) return num_changed;
    // divide by the total number of words in the corpus
    double norm = ferrum::sum(vals);
    ferrum::scalar_product(1.0/norm, &vals);
    return num_changed;
  }
  /**
   * ensure every value is >= fl
   * return the number of values changed (pre-normalization)
   */
  int floor_and_l1_norm_probs(Eigen::VectorXd& vals, double fl) {
    // make sure there are no zero counts
    int num_changed = ferrum::ensure_min(fl, vals.data(), vals.size());
    if(num_changed == 0) return num_changed;
    // divide by the total number of words in the corpus
    double norm = vals.sum();
    vals /= norm;
    return num_changed;
  }


  /**
   * Compute x = ax + by, where x, y are vectors and a, b are 
   * scalars.
   */
  template <>
  void linear_combination_in_first(std::vector< double >* const x, const std::vector< double >& y, const double a, const double b) {
    const size_t dim = x->size();
    if(dim != y.size()) {
      throw 1;
    }
    double* xp = x->data();
    cblas_dscal((const int)dim, a, xp, 1);
    const double* yp = y.data();
    cblas_daxpy((const int)dim, b, yp, 1, xp, 1);
  };

  RedirectBuffer::RedirectBuffer(std::ostream& ost) : ost_(&ost) {
    old = ost_->rdbuf(buffer.rdbuf());
  };
  RedirectBuffer::~RedirectBuffer() {
    ost_->rdbuf(old);
  }

  SmartWriter::SmartWriter() :
    console_(true), f_ptr(NULL),
    curr_file("stdout"),
    mutex_(std::shared_ptr<std::mutex>(new std::mutex)) {
  }
  SmartWriter::SmartWriter(const std::string& fname) : 
    base_(fname), console_(fname == "-"), f_ptr(NULL),
    mutex_(std::shared_ptr<std::mutex>(new std::mutex))  {
    if(console_) curr_file = "stdout";
  }

  SmartWriter::~SmartWriter() {
    if(!console_ && f_ptr != NULL) {
      f_ptr->close();
      delete f_ptr;
    }
  }

  std::mutex& SmartWriter::mutex() {
    return *(mutex_.get());
  }

  std::string SmartWriter::next_file_name(const std::string& suffix) {
    return suffix.size() > 0 ? (base_ + "." + suffix) : base_;
  }
  std::ostream& SmartWriter::get(const std::string& suffix) {
    if(!console_) {
      std::string f = next_file_name(suffix);
      if(f_ptr != NULL) {
	INFO << "Closing previously opened file " << curr_file << " and opening the new one " << f;
	f_ptr->close();
	delete f_ptr;
      }
      curr_file = f;
      f_ptr = new std::ofstream(f, std::ofstream::out);
      return *f_ptr;
    } else {
      return std::cout;
    }    
  }
  std::ostream& SmartWriter::get() {
    return get("");
  }
  std::ostream& SmartWriter::get(const int i) {
    return get(std::to_string(i));
  }
  std::string SmartWriter::base_name() {
    return base_;
  }
  std::string SmartWriter::name() {
    return curr_file;
  }

  bool SmartWriter::to_file() {
    return !console_ && this->name() != "/dev/null";
  }

  /**
   * Print various stats of the process
   */
  void print_pstats() {
    pid_t pid = getpid();
    pid_t ppid = getppid();
    INFO << "PROCESS ID: " << pid;
    INFO << "Parent Process ID: " << ppid;
    INFO << "schema++ Library built from: " << ferrum::FERRUM_GIT_SHA;
    INFO << "schema++ Library built one: " << ferrum::FERRUM_BUILD_DATE;
  }


  // template <> std::vector<double> scalar_product<double, double>(const double& x, const std::vector< double >& y);
  // template <> std::vector<double> scalar_product<int, double>(const int& x, const std::vector< double >& y);
  // template <> void scalar_product<double, double>(const double& x, std::vector< double >* y);
  // template <> void scalar_product<int, double>(const int& x, std::vector< double >* y);
}
