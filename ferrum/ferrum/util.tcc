#ifndef FERRUM_UTIL_TCC_
#define FERRUM_UTIL_TCC_

#include "ferrum/util.hpp"

#include "ferrum/cblas_cpp.hpp"
#include <mutex>
#include <thread>
#include <vector>

#include <Eigen/Dense>

namespace ferrum {
  
  /**
   * return true iff first is fully-contained in second
   */
  template <typename X, typename Y>
  inline bool contained_in(const std::pair<X, Y>& first, const std::pair<X, Y>& second) {
    return					\
      (first.first >= second.second) && 
      (first.second <= second.second) &&
      (first.first <= second.second) &&
      (first.second >= second.first);
  }

  template <typename X>
  inline std::pair<X,X> span(const std::vector<X>& v) {
    X start{};
    X end{};
    if(v.size() > 0) {
      start = v[0];
      end = v[v.size() - 1];
    }
    return std::pair<X,X>(start, end);
  }


  /**
   * Find the maximum in an M-length vector of elements.
   */
  template <typename V>
  inline typename V::value_type max(const V& counts) {
    typedef typename V::value_type T;
    T max_num{std::numeric_limits<T>::lowest()};
    for(const T& elem_count : counts) {
      if(elem_count > max_num) {
	max_num = elem_count;
      }
    }
    return max_num;
  };
  template <typename T>
  inline T max(const std::vector<T>& counts) {
    T max_num{std::numeric_limits<T>::lowest()};
    const T* const start = counts.data();
    const size_t size = counts.size();
    for(size_t i = 0; i < size; ++i) {
      T elem_count = start[i];
      if(elem_count > max_num) {
	max_num = elem_count;
      }
    }
    return max_num;
  }

  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const std::vector< T >& y) {
    const size_t dim = x->size();
    if(dim != y.size()) {
      throw 1;
    }
    for(size_t idx = 0; idx < dim; ++idx) {
      x->operator[](idx) += y[idx];
    }
  };
  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const Eigen::VectorXd& y) {
    const size_t dim = x->size();
    if(dim != y.rows()) {
      throw 1;
    }
    for(size_t idx = 0; idx < dim; ++idx) {
      x->operator[](idx) += y(idx);
    }
  }

  template <>
  void sum_in_first(std::vector< float >* const x, const std::vector< float >& y);
  template <>
  void sum_in_first(std::vector< double >* const x, const std::vector< double >& y);
  template <>
  void sum_in_first(std::vector< double >* const x, const Eigen::VectorXd& y);

  template <typename T>
  inline std::vector<T> sum(const T& x, const std::vector< T >& y) {
    const size_t dim = y.size();
    std::vector<T> res(dim, x);
    sum_in_first<T>(&res, y);
    return res;
  };

  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const std::map< int, T >& y) {
    for(const auto& p : y) {
      int idx = p.first;
      const T& prev = x->at(idx);
      x->operator[](idx) = prev + p.second;
    }
  };

  template <typename T, typename U>
  inline std::vector<T> scalar_product(const U& x, const std::vector< T >& y) {
    const size_t dim = y.size();
    std::vector<T> res(dim);
    for(size_t idx = 0; idx < dim; ++idx) {
      res[idx] = x * y[idx];
    }
    return res;
  };

  template <typename T, typename U>
  inline void scalar_product(const U& x, std::vector< T >* y) {
    const size_t dim = y->size();
    for(size_t idx = 0; idx < dim; ++idx) {
      y->operator[](idx) *= x;
    }
  };


  /**
   * Compute x = ax + by, where x, y are vectors and a, b are 
   * scalars.
   */
  template <typename T>
  inline void linear_combination_in_first(std::vector< T >* const x, const std::vector< T >& y, const T a, const T b) {
    const size_t dim = x->size();
    if(dim != y.size()) {
      throw 1;
    }
    for(size_t idx = 0; idx < dim; ++idx) {
      T val = a * (x->operator[](idx));
      x->operator[](idx) = val + (b * y[idx]);
    }
  };

  template <>
  void linear_combination_in_first(std::vector< double >* const x, const std::vector< double >& y, const double a, const double b);

  template <typename T>
  inline void atomic_incr(T* value, std::mutex* mut, const T& by) {
    std::lock_guard<std::mutex> lock(*mut);
    *value += by;
  }

}

#endif
