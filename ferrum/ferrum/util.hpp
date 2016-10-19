#ifndef FERRUM_UTIL_H_
#define FERRUM_UTIL_H_

#include "ferrum/logging.hpp"
#include <gsl/gsl_sf_log.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <map>
#include <math.h>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <type_traits> // for std::is_const
#include <unordered_map>
#include <utility> // for std::pair
#include <vector>

#include <Eigen/Dense>

#define IGNORE_RESULT(x) do { \
    (x);		      \
  } while(0)


namespace ferrum { 

  std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec);
  std::ostream& operator<<(std::ostream&& out, const std::vector<double>& vec);

  // template <typename CharT, typename Traits = std::char_traits<CharT> >
  // inline std::basic_ostream<CharT, Traits>&
  // operator<<(std::basic_ostream<CharT, Traits>&& out, const std::vector<double>& vec) {
  //   const size_t K = vec.size();
  //   out << "(";
  //   for(size_t k = 0; k < K; ++k) {
  //     out << vec.at(k);
  //     if(k + 1 < K) {
  // 	out << " ";
  //     } else {
  // 	out << ")";
  //     }
  //   }
  //   return out;
  // }

  void lower(std::string&);

  /**
   * return true iff first is fully-contained in second
   */
  template <typename X, typename Y>
  inline bool contained_in(const std::pair<X, Y>& first, const std::pair<X, Y>& second);

  template <typename X>
  inline std::pair<X,X> span(const std::vector<X>&);

  class RedirectBuffer {
  public:
    RedirectBuffer(std::ostream& ost);
    ~RedirectBuffer();
  private:
    std::stringstream buffer;
    std::streambuf* old;
    std::ostream* ost_;
  };

  class SmartWriter {
  public:
    SmartWriter();
    SmartWriter(const std::string& fname);
    virtual ~SmartWriter();
    std::ostream& get();
    std::ostream& get(const int i);
    std::ostream& get(const std::string& suffix);
    std::string base_name();
    std::string name();
    bool to_file();
    std::mutex& mutex();
  protected:
    std::string base_;
    const bool console_;
    std::ofstream* f_ptr = NULL;
    std::string curr_file;
    std::string next_file_name(const std::string& suffix);
  private:
    std::shared_ptr<std::mutex> mutex_;
  };

  /**
   * Print various stats of the process
   */
  void print_pstats();

  template <typename T>
  std::vector<T> zeros(int size_) {
    std::vector<T> zero(size_, 0);
    return zero;
  }

  template <typename T>
  inline T min(const T& a, const T& b) {
    return a < b ? a : b;
  }

  template <typename T>
  class VectorSortIndex {
  private:
    const std::vector<T>* vp_;
    int mult_;
  public:
    VectorSortIndex(const std::vector<T>* v, int mult) : vp_(v), mult_(mult) {
    }
    bool operator() (const size_t& i, const size_t& j) { 
      return (mult_ * (*vp_)[i]) < (mult_ * (*vp_)[j]);
    }
  };

  template <typename T>
  std::vector<size_t> sort_indices(const std::vector<T> &v, bool ascending) {
    const int mult = ascending ? 1 : -1;
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for(size_t i = 0; i != idx.size(); ++i) {
      idx[i] = i;
    }

    VectorSortIndex<T> vsi(&v, mult);
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(), vsi);

    return idx;
  }

  template <typename C1, typename C2>
  inline void copy(const C1* from, C2* to) {
    std::copy(from->begin(), from->end(), to->begin());
  }


  /**
   * Find the average in a vector.
   */
  template <typename T>
  inline T average(const std::vector< T >& counts) {
    T rsum = 0;
    for(const T& elem : counts) {
      rsum += elem;
    }
    return counts.size() ? (rsum / (T)counts.size()) : 0.0;
  };

  template <typename T>
  inline T median(const std::vector< T >& counts) {
    const size_t size = counts.size();
    T median{};
    if(size == 0) return median;
    std::vector<T> copy(counts);
    std::sort(copy.begin(), copy.end());
    if(size % 2 == 0) {
      median = (copy[size / 2 - 1] + copy[size / 2]) / 2;
    } else {
      median = copy[size / 2];
    }
    return median;
  };

  template <typename T>
  inline T dist(const std::vector<T>& x, const std::vector<T>& y) {
    if(x.size() != y.size()) {
      ERROR << "Sizes not equal: " << x.size() << " vs. " << y.size();
      throw 2;
    }
    double d = 0.0;
    for(size_t i = 0; i < (size_t)x.size(); ++i) {
      double diff = x[i] - y[i];
      d += (diff * diff);
    }
    return sqrt(d);
  }

  inline std::vector<double> compute_coherences(const int M,
						const std::vector< std::vector<double> >& probs,
						const std::map<int, int>& single_occur,
						const std::map< std::pair<int, int>, int >& double_occur) {
    std::vector<double> coherences;
    const int num_dists = probs.size();
    for(int dist_idx = 0; dist_idx < num_dists; ++dist_idx) {
      const std::vector<double>& curr_probs = probs[dist_idx];
      const std::vector<size_t> sorted_topic = ferrum::sort_indices(curr_probs, false);
      double run_coher = 0.0;
      for(int m = 1; m < M; ++m) { 
	const int word_at_m = sorted_topic[m];
	for(int l = 0; l < m - 1; ++l) {
	  const int word_at_l = sorted_topic[l];
	  std::pair<int,int> p( word_at_m, word_at_l );
	  double num = 1.0;
	  if(! double_occur.count(p)) {
	    p.first = word_at_l;
	    p.second = word_at_m;
	    if( !double_occur.count(p) ){
	      // do nothing
	    } else {
	      double x = double_occur.at(p);
	      num += x;
	    }
	  } else {
	    double x = double_occur.at(p);
	    num += x;
	  }
	  //double denom = (double)( single_occur.count(word_at_l) ? single_occur.at(word_at_l) : 1E-9);
	  double denom = (double)( single_occur.count(word_at_l) ? single_occur.at(word_at_l) : 0);
	  run_coher += gsl_sf_log( num / denom );
	}
      }
      coherences.push_back(run_coher);
    }
    return coherences;
  }

  template <typename T>
  inline std::vector<T> sum(const std::vector< T >& x, const std::vector< T >& y) {
    if(x.size() != y.size()) {
      throw 1;
    }
    const size_t dim = x.size();
    std::vector<T> res(dim);
    for(size_t idx = 0; idx < dim; ++idx) {
      res[idx] = x[idx] + y[idx];
    }
    return res;
  };

  template <typename T>
  inline T dot(const std::vector< T >& x, const std::vector< T >& y) {
    if(x.size() != y.size()) {
      ERROR << "x.size " << x.size() << " != y.size " << y.size();
      throw 1;
    }
    const size_t dim = x.size();
    T res{};
    for(size_t idx = 0; idx < dim; ++idx) {
      T prod = x[idx] * y[idx]; 
      res = res + prod;
    }
    return res;
  };

  template <typename T>
  inline void sum(const T& x, std::vector< T >* y) {
    const size_t dim = y->size();
    for(size_t idx = 0; idx < dim; ++idx) {
      (*y)[idx] += x;
    }
  };

  /**
   * Compute x = ax + by, where x, y are vectors and a, b are 
   * scalars.
   */
  template <typename T>
  inline void linear_combination_in_first(std::vector< T >* const x, const std::vector< T >& y, const T a, const T b);

  /**
   * Compute x = x + y, where x, y are vectors.
   */
  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const std::vector< T >& y);
  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const Eigen::VectorXd& y);
  template <typename T>
  inline void sum_in_first(std::vector< T >* const x, const std::map< int, T >& y);

  /**
   * Compute x = a + x, where x is a vector and a is a scalar.
   */
  template <typename T>
  inline std::vector<T> sum(const T& x, const std::vector< T >& y);

  template <typename T>
  inline int ensure_min(const T& val, std::vector< T >* const x) {
    const size_t dim = x->size();
    int num_changed = 0;
    for(size_t idx = 0; idx < dim; ++idx) {
      const T& xval = (*x)[idx];
      if(xval < val) {
	(*x)[idx] = val;
	num_changed++;
      }
    }
    return num_changed;
  };

  template <typename T>
  inline int ensure_min(const T& val, T* const x, size_t dim) {
    int num_changed = 0;
    for(size_t idx = 0; idx < dim; ++idx) {
      const T& xval = x[idx];
      if(xval < val) {
	x[idx] = val;
	num_changed++;
      }
    }
    return num_changed;
  };

  template <typename T>
  inline std::vector<T> column(const std::vector<std::vector<T> >& matrix, const size_t& col_idx) {
    std::vector<T> slice;
    for(const std::vector<T>& row : matrix) {
      slice.push_back(row[col_idx]);
    }
    return slice;
  }

  template <typename T>
  inline bool vectors_equal(const std::vector<T>& x, const std::vector<T>& y) {
    bool eq = true;
    if(x.size() != y.size()) {
      ERROR << "comparing containers of different sizes";
      return false;
    }
    const size_t size = x.size();
    for(size_t i = 0; i < size; ++i) {
      eq &= (x[i] == y[i]);
      if(!eq) {
	return false;
      }
    }
    return true;
  }

  template <typename M >
  inline bool maps_equal(const M& x, const M& y) {
    bool eq = true;
    if(x.size() != y.size()) {
      ERROR << "comparing containers of different sizes";
      return false;
    }
    for(typename M::const_iterator it = x.begin();
	it != x.end();
	++it) {
      typename M::const_iterator yit = y.find( it->first );
      eq &= ( yit != y.end() );
      if(!eq) {
	return false;
      }
      eq &= ( yit->second == it->second);
      if(!eq) {
	return false;
      }
    }
    return true;
  }

  template <typename T>
  inline void product(const T& val, std::vector< T >* const x) {
    const int dim = x->size();
    for(int idx = 0; idx < dim; ++idx) {
      (*x)[idx] *= val;
    }
  };

  template <typename T>
  inline void invert(std::vector< T >* const x) {
    const int dim = x->size();
    for(int idx = 0; idx < dim; ++idx) {
      const T& val = x->operator[](idx);
      x->operator[](idx) = 1.0/val;
    }
  };

  template <typename V>
  inline void zero_out(V* const x) {
    typename V::iterator it = x->begin();
    for(;it != x->end(); ++it) {
      *it = 0.0;
    }
  };

  // template <typename V>
  // inline void exp(V* const x) {
  //   //typedef typename V::value_type T;
  //   typename V::iterator it = x->begin();
  //   // if(std::is_const< decltype(it) >::value) {
  //   // 	for(;it != x->end(); ++it) {
  //   // 	  T temp = *it;
  //   // 	  x->erase(it);
  //   // 	  x->insert(it, gsl_sf_exp(temp));
  //   // 	}
  //   // } else {
  //   for(;it != x->end(); ++it) {
  // 	*it = mathops::exp(*it);
  //   }
  //   //      }
  // };
  // template <typename V>
  // inline void log(V* const x) {
  //   for(typename V::iterator it = x->begin();
  // 	  it != x->end(); ++it) {
  // 	*it = gsl_sf_log(*it);
  //   }
  // };
  // template <typename V>
  // inline V exp(const V& x) {
  //   V res(x);
  //   exp(&res);
  //   return res;
  // };
  // template <typename V>
  // inline V log(const V& x) {
  //   V res(x);
  //   log(&res);
  //   return res;
  // };

  double variance(const Eigen::VectorXd&);

  template <typename V>
  inline void root_sq(V* const x) {
    typedef typename V::value_type T;
    typename V::iterator it = x->begin();
    for(;it != x->end(); ++it) {
      T val = *it;
      *it = sqrt(val);
    }
  };

  template <typename V>
  inline void square(V* const x) {
    typedef typename V::value_type T;
    typename V::iterator it = x->begin();
    for(;it != x->end(); ++it) {
      T val = *it;
      *it = val * val;
    }
  };
  template <typename V>
  inline V square(const V& x) {
    V res(x);
    square(&res);
    return res;
  };
  template <typename V>
  inline void cube(V* const x) {
    typedef typename V::value_type T;
    typename V::iterator it = x->begin();
    for(;it != x->end(); ++it) {
      T val = *it;
      *it = val * val * val;
    }
  };
  template <typename V>
  inline V cube(const V& x) {
    V res(x);
    cube(&res);
    return res;
  };
  template <typename V>
  inline void quartic(V* const x) {
    typedef typename V::value_type T;
    typename V::iterator it = x->begin();
    for(;it != x->end(); ++it) {
      const T val = *it;
      *it = val * val * val * val;
    }
  };
  template <typename V>
  inline V quartic(const V& x) {
    V res(x);
    quartic(&res);
    return res;
  };

  template <typename U>
  std::vector<double> scalar_product(const U& x, const std::vector< double >& y);
  template <typename T, typename U>
  std::vector<T> scalar_product(const U& x, const std::vector< T >& y);

  template <typename U>
  void scalar_product(const U& x, std::vector< double >* y);
  template <typename T, typename U>
  void scalar_product(const U& x, std::vector< T >* y);

  void product(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&,
	       const std::vector<double>&,
	       std::vector<double>*);
  void product(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&,
	       const double*,
	       double*);

  template <typename T>
  inline void product(const std::vector< std::vector< T > >& mat,
		      const std::vector< T > &vec,
		      std::vector< T >* acc) {
    const size_t c = mat[0].size();
    const size_t vr = vec.size();
    const size_t mr = mat.size();
    const size_t ar = acc->size();
    if(vr != c || mr != ar) {
      ERROR << "mismatched dimension of (" << mr << " x " << c << ") * (" << vr << " x 1) == (" << ar << " x 1)";
      throw 5;
    }
    for(size_t i = 0; i < mr; ++i) {
      for(size_t j = 0; j < c; ++j) {
	acc->operator[](i) += mat[i][j]*vec[j];
      }
    }
  }
  template <typename T>
  inline void transpose_product(const std::vector< std::vector< T > >& mat,
				const std::vector< T > &vec,
				std::vector< T >* acc) {
    const size_t c = mat[0].size();
    const size_t vr = vec.size();
    const size_t mr = mat.size();
    const size_t ar = acc->size();
    if(ar != c || mr != vr) {
      ERROR << "mismatched dimension of (1 x " << vr << ") * (" << mr << " x " << c << ") == (" << ar << " x 1)";
      throw 5;
    }
    for(size_t i = 0; i < mr; ++i) {
      const T& vi = vec[i];
      for(size_t j = 0; j < c; ++j) {
	acc->operator[](j) += mat[i][j]*vi;
      }
    }
  }

  // template <typename T>
  // inline void transpose(std::vector< std::vector< T > >* mat);

  template <typename T>
  inline void set(std::vector<T>* container, const T& value) {
    std::fill(container->begin(), container->end(), value);
  }

  template <typename T>
  inline void atomic_incr(T* value, std::mutex* mut, const T& by);

  

  /**
   * Find the maximum in an M-length vector of elements.
   */
  template <typename V> typename V::value_type max(const V& counts);
  /**
   * Find the maximum in an M-length vector of elements.
   */
  template <typename V> V max(const std::vector<V>& counts);

  template <typename T>
  inline size_t arg_max(const std::vector<T>& counts) {
    T max_num{std::numeric_limits<T>::lowest()};
    size_t which = 0, i = 0;
    for(const T& elem_count : counts) {
      if(elem_count > max_num) {
	max_num = elem_count;
	which = i;
      }
      i++;
    }
    return which;
  };

  /**
   * Find the maximum in an MxN matrix of elements.
   */
  template <typename T>
  inline T max(const std::vector< std::vector< T > >& counts) {
    T max_num{std::numeric_limits<T>::lowest()};
    for(const std::vector<T>& row_counts : counts) {
      T r_max = max(row_counts);
      if(r_max > max_num) {
	max_num = r_max;
      }
    }
    return max_num;
  };

  /**
   * Find the maximum in an MxN matrix of elements.
   */
  template <typename T>
  inline T max(const std::vector<std::vector< std::vector< T > > >& counts) {
    T max_num{std::numeric_limits<T>::lowest()};
    for(const std::vector<std::vector<T> >& row_counts : counts) {
      T r_max = max<T>(row_counts);
      if(r_max > max_num) {
	max_num = r_max;
      }
    }
    return max_num;
  };

  /**
   * Find the maximum row elements in an MxN matrix of elements.
   */
  template <typename T>
  inline std::vector<T> row_max(const std::vector< std::vector< T > >& counts) {
    const size_t M = counts.size();
    std::vector<T> maxes(M);
    size_t m = 0;
    for(const std::vector<T>& row_counts : counts) {
      maxes[m++] = max(row_counts);
    }
    return maxes;
  };
  template <typename T>
  inline std::vector<T> row_arg_max(const std::vector< std::vector< T > >& counts) {
    const size_t M = counts.size();
    std::vector<T> maxes(M);
    size_t m = 0;
    for(const std::vector<T>& row_counts : counts) {
      maxes[m++] = arg_max(row_counts);
    }
    return maxes;
  };

  /**
   * Find the maximum in an MxN matrix of elements.
   */
  template <typename T>
  inline std::vector<T> column_max(const std::vector< std::vector< T > >& counts) {
    const int K = counts[0].size();
    std::vector<T> maxes(K);
    for(const std::vector<T>& row_counts : counts) {
      for(int k = 0; k < K; ++k) {
	if(row_counts[k] > maxes[k]) {
	  maxes[k] = row_counts[k];
	}
      }
    }
    return maxes;
  };

  /**
   * Find the maximum in an MxN matrix of elements.
   */
  template <typename T>
  inline std::vector<T> column_max(const std::vector<std::vector< std::vector< T > > >& counts) {
    const int K = counts[0][0].size();
    std::vector<T> maxes(K);
    for(const std::vector<std::vector<T> >& row_counts : counts) {
      std::vector<T> m = column_max<T>(row_counts);
      for(int k = 0; k < K; ++k) {
	if(m[k] > maxes[k]) {
	  maxes[k] = m[k];
	}
      }
    }
    return maxes;
  };

  /**
   * Given a count matrix of D x K, return a histogram
   * recording, for each k <= K, the number of times k appeared.
   */
  inline std::vector<std::vector<int> > histogram(const std::vector<std::vector< int > >& counts) {
    const int K = counts[0].size();
    const int num_rows = counts.size();
    const int matrix_max = ferrum::max(counts);
    std::vector< std::vector< int > > hist(K, std::vector<int>(matrix_max + 1, 0));
    for(int k = 0; k < K; ++k) {
      for(int i = 0; i < num_rows; ++i) {
	++hist[k][ counts[i][k] ];
      }
    }
    return hist;
  };

  /**
   * Given a count matrix of D x K, return a histogram
   * recording, for each k <= K, the number of times k appeared.
   */
  inline std::vector<std::vector<int> > mf_histogram(const std::vector< std::vector<std::vector< int > > >& counts) {
    const int K = counts[0][0].size();
    const int num_rows = counts.size();
    const int num_cols = counts[0].size();
    const int matrix_max = ferrum::max(counts);
    std::vector< std::vector< int > > hist(K, std::vector<int>(matrix_max + 1, 0));
    for(int k = 0; k < K; ++k) {
      for(int i = 0; i < num_rows; ++i) {
	for(int j = 0; j < num_cols; ++j) {
	  ++hist[k][ counts[i][j][k] ];
	}
      }
    }
    return hist;
  };

  template <typename V> 
  inline typename V::value_type sum(const V& counts) {
    typedef typename V::value_type T;
    T sum = {};
    for(const T& x : counts) {
      sum += x;
    }
    return sum;
  }

  inline std::map<int,int> sparse_histogram(const std::vector< int >& counts) {
    std::map<int, int> hist_;
    for(const int x : counts) {
      ++hist_[x];
    }
    return hist_;
  }

  inline std::vector<int> marginals(const std::vector<std::vector< int > >& counts) {
    const int num_rows = counts.size();
    std::vector<int> marginals_(num_rows, 0);
    int i = 0;
    for(const std::vector<int>& row_counts : counts) {
      marginals_[i++] = ferrum::sum(row_counts);
    }
    return marginals_;
  }

  inline std::vector<int> marginals(const std::vector<std::vector<std::vector< int > > >& counts) {
    std::vector<int> marginals_;
    for(const std::vector<std::vector<int> >& row_counts : counts) {
      for(const std::vector<int>& col_counts : row_counts) {
	marginals_.push_back(ferrum::sum(col_counts));
      }
    }
    return marginals_;
  }

  /**
   * 
   */
  inline std::vector<int> marginal_histogram(const std::vector<std::vector< int > >& counts) {
    const std::vector<int> marginals_ = ferrum::marginals(counts);
    const int max_marginal = ferrum::max(marginals_);
    std::vector<int> lengths(max_marginal + 1, 0);
    for(const int marg : marginals_) {
      ++lengths[marg];
    }
    return lengths;
  };

  /**
   * ensure every value is >= fl
   * return the number of values changed (pre-normalization)
   */
  int floor_and_l1_norm_probs(std::vector<double>& vals, double fl);
  /**
   * ensure every value is >= fl
   * return the number of values changed (pre-normalization)
   */
  int floor_and_l1_norm_probs(Eigen::VectorXd& vals, double fl);
  /**
   * Find the Frobenius norm of a (dense) matrix. 
   * The matrix cannot be ragged.
   */
  template <typename T>
  inline double frobenius_norm(const std::vector< std::vector< T > >& matrix) {
    double fn = 0.0;
    for(const auto& row : matrix) {
      for(const auto& elem : row) {
	fn += (elem * elem);
      }
    }
    return sqrt(fn);
  };

  template <typename T>
  inline double l2_norm(const std::vector<T>& vec) {
    double fn = 0.0;
    for(const auto& elem : vec) {
      fn += (elem*elem);
    }
    return sqrt(fn);
  }


  template <typename T>
  class pair_hash { 
  public:
    const size_t operator()(const std::pair<T, T>& pair_) const {
      return std::hash<T>()(pair_.first) ^ std::hash<T>()(pair_.second);
    }
  };

  template <typename K, typename V > using pair_map =
    std::unordered_map< const std::pair< K, K >, V , ferrum::pair_hash<K> >;
  template <typename K> using pair_icount =
    std::unordered_map< const std::pair< K, K >, int, ferrum::pair_hash<K> >;

  inline void print_1d(const Eigen::VectorXd& vec) {
    for(int i = 0; i < vec.rows(); ++i) {
      std::cout << std::setprecision(8) << vec(i) << " ";
      //	printf("%.8f ", elem);
    }
    //printf("\n");
    std::cout << std::endl;
  }

  template <typename V>
  inline void print_1d(const V& vec) {
    for(const auto& elem: vec) {
      std::cout << std::setprecision(8) << elem << " ";
      //	printf("%.8f ", elem);
    }
    //printf("\n");
    std::cout << std::endl;
  }
  template <typename T>
  inline void print_2d(const std::vector<std::vector<T> >& vec) {
    for(const std::vector<T>& rvec : vec) {
      for(const T& elem: rvec) {
	std::cout << elem << ' ';
      }
      std::cout << std::endl;
    }
  }
  
  template <typename T>
  inline void print_2d_distribution(const std::vector<std::vector<T> >& vec, std::ostream & outter = std::cout) {
    for(const std::vector<T>& rvec : vec) {
      std::stringstream stream;
      for(const T& elem: rvec) {
	stream << elem << ' ';
      }
      outter << stream.str();
      outter << std::endl;
    }
  }

  template <typename F>
  inline F str_num_value(const std::string& str) {
    F res;
    std::istringstream convert(str);
    if( ! (convert >> res)) {
      ERROR << "Could not parse string " << str << " into a numeric type";
      res = (F)0.0;
    }
    return res;
  }
  template <typename F>
  inline F str_to_value(const std::string& str) {
    F res;
    std::istringstream convert(str);
    if( ! (convert >> res)) {
      ERROR << "Could not parse string " << str << " into a numeric type";
      throw 2;
    }
    return res;
  }

  template <typename Container>
  inline void to_stringstream(const Container& container, std::stringstream& ss) {
    size_t i = 0;
    const size_t size = container.size();
    for(const auto& item : container) {
      ss << item;
      ++i;
      if(i < size) {
	ss << " ";
      }
    }
  }

  template <typename T>
  inline void check_positive(const std::vector< T >& vec) {
    int which = 0;
    for(const T& elem : vec) {
      if(elem <= 0) {
	ERROR << "Element[" << which << "] = " << elem << " is not positive";
      }
      which++;
    }
  }

  template <typename T>
  inline void check_positive(const std::vector<std::vector<T> >& vec) {
    for(const std::vector<T>& rvec : vec) {
      check_positive<T>(rvec);
    }
  }

  template <typename T>
  inline void check_positive(const std::vector<std::vector<std::vector<T> > >& vec) {
    for(const std::vector< std::vector<T > >& rvec : vec) {
      check_positive<T>(rvec);
    }
  }

  inline void* MALLOC(size_t size) {
    void* ptr = malloc(size);
    if(ptr == NULL) {
      throw 5;
    }
    return ptr;
  }

  template <class T> inline void fill_array1d(T* array, T value,
					      const int size) {
    DEBUG << "filling " << array << " with " << size << " copies of " << value;
    for(int i = 0; i < size; i++) {
      array[i] = value;
    }
  }

  template <typename X> inline void allocate_1d(X* &array, const int size) {
    DEBUG << "making 1D array @ " << &array << " of size " << size;
    array = (X*)MALLOC(sizeof(X) * size);
  }

  template <typename X> inline void allocate_2d(X** &array, const int num_rows, const int* num_cols) {
    DEBUG << "making 2D array @ " << &array << " of row size " << num_rows;
    //array = (X**)MALLOC(sizeof(X*) * num_rows);
    allocate_1d<X*>(array, num_rows);
    for(int i = 0; i < num_rows; i++) {
      DEBUG << "\trow " << i << ", num cols: " << num_cols[i];
      DEBUG << " @ " << array[i];
      DEBUG << *(array+i);
      allocate_1d<X>(array[i], num_cols[i]);
    }
    DEBUG << "done";
  }

  template <typename X> inline const X** transfer_2d(const X** into,
						     X** source) {
    const int r_len = sizeof(source)/sizeof(X*);
    into = (const X**)malloc(sizeof(X*) *  r_len);
    for(int r = 0; r < r_len; r++) {
      const int c_len = sizeof(into[r])/sizeof(X);
      into[r] = (const X*)malloc(sizeof(X) * c_len);
      for(int c = 0; c < c_len; c++) {
	into[r][c] = source[r][c];
      }
    }
    return into;
  }
  template <typename X> inline const X*& transfer_1d(X* source) {
    X* into;
    const int r_len = sizeof(source)/sizeof(X);
    DEBUG << "For transfer: allocating array @ " << &into << " of size " << r_len;
    allocate_1d(into, r_len);
    for(int r = 0; r < r_len; r++) {
      into[r] = source[r];
    }
    return into;
  }

  template <typename X> inline X* make_1d(const int size, const X value) {
    X* into;
    allocate_1d(into, size);
    fill_array1d<X>(into, value, size);
    return into;
  }


  // Inline functions

  template <typename X> inline X& get_1d_ref(X arr) {
    return &arr;
  }

  /////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////

  template <typename T> class Vector {
  protected:
    std::vector< T > array;
  public:
    Vector<T>() {
    }
    T& operator[](std::size_t idx) {
      return array[idx];
    }
    void add_row() {
      array.push_back( T() );
    }
    void add_row(const int size, T val) {
      const int curr_row = array.size();
      array[curr_row].push_back(val);
    }
  };

  template <typename T> 
  class Vector2D {
  private:
    std::vector< std::vector<T> > array;
  public:
    Vector2D<T>() {
      DEBUG << "calling vector2d constructor";
    }
    
    // virtual ~Vector2D() {
    //   DEBUG << "calling vector2d destructor" << std::endl;
    // }

    std::vector<T>& operator[](std::size_t idx) {
      return array[idx];
    }

    void add_row() {
      std::vector<T>* vec = new std::vector<T>();
      array.push_back( *vec );
    }
    void add_row(const int size, T* arr) {
      const int curr_row = array.size();
      DEBUG << "Curr_row = " << curr_row;
      add_row();
      for(int i = 0; i < size; i++) {
	array[curr_row].push_back(arr[i]);
      }
    }
  };
}

#define CAPTURE_STDOUT ferrum::RedirectBuffer ____stdout_capturer__with_very_unique_name(std::cout)
#define CAPTURE_STDERR ferrum::RedirectBuffer ____stderr_capturer__with_very_unique_name(std::cerr)

#endif

#include "ferrum/util.tcc"
