#include "ferrum/logging.hpp"
#include "ferrum/optimize.hpp"
#include <math.h>
#include <cstring>
#include <memory>

namespace optimize {
  std::shared_ptr<Minimizer>
  Minimizer::make(OptimizationMethod opt_how, size_t support) {
    std::shared_ptr<Minimizer> min;
    if(opt_how == optimize::OptimizationMethod::LBFGS) {
#ifndef PUBLIC_UPF
      min =
	std::shared_ptr<optimize::LibLBFGSMinimizer>
	(new optimize::LibLBFGSMinimizer(support));
#else
      throw 10;
#endif
    } else if(opt_how == optimize::OptimizationMethod::ADAGRAD) {
      min =
	std::shared_ptr<optimize::AdaGrad>(new optimize::AdaGrad(support));
    } else {
      ERROR << "Unknown optimization method " << opt_how;
    }
    return min;
  }

  GSLMinimizer::GSLMinimizer(const int num_dim) : minimizer_type_(gsl_multimin_fdfminimizer_vector_bfgs2), minimizer_(gsl_multimin_fdfminimizer_alloc(minimizer_type_, num_dim)), gradient_tolerance_(1E-3), tolerance_(1E-4), step_size_(1E-2), /*num_dim_(num_dim),*/ num_steps_(100) {
  }
  GSLMinimizer::~GSLMinimizer() {
    gsl_multimin_fdfminimizer_free(minimizer_);
  }
  void GSLMinimizer::gradient_tolerance(double gt) {
    gradient_tolerance_ = gt;
  }
  void GSLMinimizer::tolerance(double tol) {
    tolerance_ = tol;
  }
  void GSLMinimizer::step_size(double ss) {
    step_size_ = ss;
  }
  void GSLMinimizer::num_steps(int num) {
    num_steps_ = num;
  }
  double GSLMinimizer::gradient_tolerance() {
    return gradient_tolerance_;
  }
  double GSLMinimizer::tolerance() {
    return tolerance_;
  }
  double GSLMinimizer::step_size() {
    return step_size_;
  }
  int GSLMinimizer::num_steps() {
    return num_steps_;
  }

  double GSLMinimizer::value() {
    return minimizer_->f;
  }

  GSLVector::GSLVector(const gsl_vector* gvec) : num_dim_(gvec->size), vec_(gsl_vector_alloc(num_dim_)) {
    gsl_vector_memcpy(vec_, const_cast<gsl_vector*>(gvec));
  }
  GSLVector::GSLVector(const int size) : num_dim_(size), vec_(gsl_vector_alloc(num_dim_)) {
  }
  GSLVector::~GSLVector() {
    gsl_vector_free(vec_);
  }
  std::ostream& operator<< (std::ostream& stream, const GSLVector& vec) {
    stream << "gsl_vector[";
    for(size_t i = 0; i < (size_t)vec.num_dim_; ++i) {
      stream << gsl_vector_get(vec.vec_, i);
      if(i+1 < (size_t)vec.num_dim_) {
	stream << ", ";
      }
    }
    stream << "]";
    return stream;
  }
  void GSLVector::update(gsl_vector* update) {
    gsl_vector_memcpy(vec_, const_cast<gsl_vector*>(update));
  }
  gsl_vector* GSLVector::get() {
    return vec_;
  }

  double GSLVector::dist(gsl_vector* x0, gsl_vector* x1) {
    double d = 0.0;
    for(size_t i = 0; i < (size_t)x0->size; ++i) {
      double diff = gsl_vector_get(x1, i) - gsl_vector_get(x0, i);
      d += (diff * diff);
    }
    return sqrt(d);
  }

#ifndef PUBLIC_UPF
  void LibLBFGSVector::init_vec() {
    vec_ = lbfgs_malloc(num_dim_);
    if(vec_ == NULL) {
      ERROR << "LibLBFGSVector cannot alloc enough memory (" << num_dim_ << ")";
      throw 4;
    }
  }

  LibLBFGSVector::LibLBFGSVector(const int size) : num_dim_(size) {
    init_vec();
  }
  LibLBFGSVector::~LibLBFGSVector() {
    if(vec_ != NULL) {
      lbfgs_free(vec_);
    }
  }
  lbfgsfloatval_t* LibLBFGSVector::get() {
    return vec_;
  }
#endif
  Minimizer::~Minimizer() {
  }

#ifndef PUBLIC_UPF
  LibLBFGSMinimizer::~LibLBFGSMinimizer() {
  }

  LibLBFGSMinimizer::LibLBFGSMinimizer(const int num_dim) : num_dim_(num_dim) {
    lbfgs_parameter_init(&params_);
    params_.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
  }
  LibLBFGSMinimizer::LibLBFGSMinimizer() {
    lbfgs_parameter_init(&params_);
    params_.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
  }

  int
  LibLBFGSMinimizer::minimize(Function* func, std::vector<double>& point) {
    lbfgsfloatval_t func_val;
    LibLBFGSVector vec(point);
    lbfgsfloatval_t* raw_lbfgs_point = vec.get();
    LibLBFGSFunction* func_ptr = static_cast<LibLBFGSFunction*>(func);
    int status = lbfgs(num_dim_, raw_lbfgs_point, &func_val,
		       func_ptr->eval, func_ptr->progress,
		       func_ptr->params, &params_);
    vec.to_container(point);
    return status;
  }
#endif
  AdaGrad::AdaGrad(size_t support) :
    AdaGrad(support, 0.1, 1E-8) {
  }
  AdaGrad::AdaGrad(size_t support, double eta, double epsilon) :
    support_(support),
    sum_sq_grad_(std::vector<double>(support_, 0.0)),
    eta_(eta),
    epsilon_(epsilon),
    max_steps_(1),
    f_eval_diff_(1E-6) {
    if(epsilon_ < 1E-12) {
      WARN << "Epsilon (" << epsilon_ << ") is < 1E-12: numerical instability may occur";
    }
    if(epsilon_ < 0) {
      ERROR << "Epsilon cannot be < 0 (currently " << epsilon_ << ")";
      throw 20;
    }
    if(eta_ < 1E-12) {
      WARN << "eta (" << eta_ << ") is < 1E-12: numerical instability may occur";
    }
    if(eta_ < 0) {
      ERROR << "eta should not be < 0 (currently " << eta_ << "); this implementation of AdaGrad *minimizes*";
    }
  }
  AdaGrad::~AdaGrad() {
  }
  int AdaGrad::minimize(Function* func,
			std::vector<double>& point) {
    int subres = minimize(func, point, (unsigned int)max_steps_);
    return subres;
  }
  double AdaGrad::eta() const {
    return eta_;
  }
  double AdaGrad::epsilon() const {
    return epsilon_;
  }
  const std::vector<double>&
  AdaGrad::sgs() const {
    return sum_sq_grad_;
  }
  void AdaGrad::rollback(const std::vector<double>& grad,
			std::vector<double>& point) {
    double* ss = sum_sq_grad_.data();
    double* p = point.data();
    const double* g = grad.data();
    for(size_t i = 0; i < support_; ++i) {
      double scale = epsilon_ + ss[i];
      scale = sqrt(scale);
      scale = eta_/scale;
      p[i] += scale*g[i];
    }
  }
  void AdaGrad::grad_op(const std::vector<double>& grad,
			std::vector<double>& point) {
    double* ss = sum_sq_grad_.data();
    double* p = point.data();
    const double* g = grad.data();
    for(size_t i = 0; i < support_; ++i) {
      const double gi = g[i];
      ss[i] += gi*gi;
      double scale = epsilon_ + ss[i];
      scale = sqrt(scale);
      scale = eta_/scale;
      p[i] -= scale*gi;
    }
    // std::vector<double> grad = ferrum::square(in_grad);
    // //ferrum::sum_in_first(&sum_sq_grad_, grad);
    // cblas_daxpy(support_, 1.0, grad.data(), 1,
    // 		sum_sq_grad_.data(), 1);
    // cblas_dcopy(support_, sum_sq_grad_.data(), 1,
    // 		grad.data(), 1);
    // ferrum::sum(epsilon_, &grad);
    // ferrum::root_sq(&grad);
    // ferrum::invert(&grad);
    // cblas_dscal(support_, -1.0*eta_, grad.data(), 1);
    
    // //ferrum::linear_combination_in_first(&point, grad, 1.0, -1.0*eta_);
    // catlas_daxpby(support_,
    // 		  -1.0*eta_, grad.data(), 1,
    // 		  1.0, point.data(), 1);
  }
  int AdaGrad::minimize(Function* func,
			std::vector<double>& point,
			unsigned int steps) {
    bool bad = false;
    AdaGradFunction* func_ptr = static_cast<AdaGradFunction*>(func);
    if(func_ptr == NULL) {
      ERROR << "Function wrapper pointer is null";
      bad = true;
    } else if(func_ptr->fdf == NULL) {
      ERROR << "Function & gradient pointer is null";
      bad = true;
    } else if(func_ptr->f == NULL) {
      ERROR << "Function pointer is null";
      bad = true;
    }
    if(support_ == 0) {
      ERROR << "Calling AdaGrad::minimize with support == 0";
      bad = true;
    } else if(point.size() != support_) {
      ERROR << "Calling AdaGrad::minimize with support " << support_ << " != point size " << point.size();
      bad = true;
    }
    if(bad) throw 10;
    bool converged = false;
    double prev_ll = std::numeric_limits<double>::max();
    unsigned int step_number = 0;
    optimize::Minimizer::Result res = optimize::Minimizer::Result::SUCCESS;
    while(! converged ) {
      std::vector<double> grad(support_, 0.0);
      prev_ll = func_ptr->f(func_ptr->params, point.data(), support_);
      func_ptr->fdf(func_ptr->params, point.data(), grad.data(), support_);
      grad_op(grad, point);
      double curr_ll = func_ptr->f(func_ptr->params, point.data(), support_);
      if((++step_number) >= steps) {
	res = optimize::Minimizer::Result::MAX_ITERATIONS;
	//INFO << "HERE, " << res;
	converged = true;
      }
      const bool f_term = std::abs(curr_ll - prev_ll) < f_eval_diff_;
      converged = converged || f_term;
      if(f_term) {
	double frob = 0.0;
	for(double x : grad) {
	  frob += x*x;
	}
	frob = sqrt(frob);
	if(frob > f_eval_diff_) {
	  res = optimize::Minimizer::Result::F_EVAL_LIMIT;
	}
      }
      if(curr_ll > prev_ll) {
	// rollback
	rollback(grad, point);
	res = optimize::Minimizer::Result::OVERSHOT;
      }
      //INFO << prev_ll << " --> " << curr_ll << "(" << std::abs(curr_ll - prev_ll) << ")";
    }
    return res;
  }

  AdaGradFunction::AdaGradFunction(AdaGrad::adagrad_feval_t f_eval,
				   AdaGrad::adagrad_fgrad_t fgrad_eval,
				   void* fparams) :
    f(f_eval),
    fdf(fgrad_eval),
    params(fparams) {
  }
  double AdaGradFunction::value(const std::vector<double>& point) {
    double x = f(params, point.data(), (int)(point.size()));
    return x;
  }
  void AdaGradFunction::grad(const std::vector<double>& point, void* result) {
    fdf(params, point.data(), (double*)result, (int)(point.size()));
  }
#ifndef PUBLIC_UPF
  LibLBFGSFunction::LibLBFGSFunction() :
    LibLBFGSFunction(NULL, NULL, NULL) {
  }
  LibLBFGSFunction::LibLBFGSFunction(lbfgs_evaluate_t fn_g_eval,
				     lbfgs_progress_t prog,
				     void* fparams) :
    eval(fn_g_eval),
    progress(prog),
    params(fparams) {
  }
  double LibLBFGSFunction::value(const std::vector<double>& point) {
    typedef lbfgsfloatval_t flt;
    std::vector<double> g(point.size(), 0.0);
    double x = eval(params, (const flt*)point.data(),
		    (flt*)g.data(), (int)(point.size()), 1.0);
    return x;
  }
  void LibLBFGSFunction::grad(const std::vector<double>& point, void* result) {
    typedef lbfgsfloatval_t flt;
    eval(params, (const flt*)point.data(),
	 (flt*)result, (int)(point.size()), 1.0);
  }
  double LibLBFGSNoOp::value(const std::vector<double>& point) {
    return 0.0;
  }
  void LibLBFGSNoOp::grad(const std::vector<double>& point, void* result) {
    double* g = (double*)result;
    for(size_t k = 0; k < point.size(); ++k) {
      *(g+k) = 0.0;
    }
  }
#endif

  double Function::linesearch
  (
   const std::vector<double>& point,
   double control_armijo,
   double control_wolfe,
   double init_step,
   Function::Want want,
   bool warn
   ) {
    // compute the initial value of the objective
    const double init_val = value(point);
    const size_t dim = point.size();
    // compute the initial gradient value
    std::vector<double> grad_val(dim, 0.0);
    this->grad(point, (void*)grad_val.data());
    // compute the initial gradient direction
    std::vector<double> g_direction(grad_val);
    ferrum::scalar_product(1.0/ferrum::l2_norm(grad_val), &g_direction);
    double dir_correction = (want == Function::Want::LOWER) ? -1.0 : 1.0;
    double init_grad_curve = dir_correction * ferrum::dot(g_direction, grad_val);
    /////////////////////////////////////////////////
    double curr_val =
      (want == Function::Want::LOWER) ?
      std::numeric_limits<double>::max() :
      std::numeric_limits<double>::lowest();
    double factors[] = {2.1, 0.5};
    // indexes into factors[]
    //int direction = (want == Function::Want::LOWER) ? 1 : 0;
    int direction = 0;
    double tolerances[] = {1E-8, 1E3};
    double step = init_step;
    double max_steps = 40;
    std::vector<double> trial_point(dim, 0.0);
    int iter = 0;
    double best_val = init_val;
    double argbest_step = 0.0;
    while(true) {
      // get the perturbed point: x = x_0 + step * gradient_direction
      {
	std::copy(point.begin(), point.end(), trial_point.begin());
	ferrum::linear_combination_in_first(&trial_point, g_direction,
					    1.0, dir_correction * step);
      }
      // get the function
      curr_val = value(trial_point);
      this->grad(trial_point, (void*)grad_val.data());
      ++iter;

      DEBUG << "Initial value: " << init_val << "; current value: " << curr_val << (direction == 0 ? "; increasing" : "; decreasing") << "; iteration " << iter;
#ifdef FERRUM_DEBUG
      ferrum::print_1d(trial_point);
#endif

      // if(curr_val > init_val) {
      // 	// otherwise, we *are* making good progress
      // 	if(curr_val > best_val) {
      // 	  argbest_step = step;
      // 	  best_val = curr_val;
      // 	}
      // 	if(iter == max_steps) {
      // 	  break;
      // 	}
      // 	goto increase_step;
      // }			       
      
      if(want == Function::Want::LOWER) {
	if(curr_val < best_val) {
      	  argbest_step = step;
      	  best_val = curr_val;
      	}
	if(curr_val > init_val + (control_armijo * step * init_grad_curve)) {
	  direction = 1;
	  TRACE << "insufficient progress; decreasing now";
	} else {
	  // if we've used up 3/4 of our iterations, then we'll take just satisfying the
	  // Armijo conditions
	  if(iter >= (int)(0.75 * max_steps)) {
	    break;
	  }
	  // now, what about the Wolfe (curvature) conditions?
	  double curr_grad_curv =
	    dir_correction * ferrum::dot(g_direction, grad_val);
	  if(curr_grad_curv < control_wolfe * init_grad_curve) {
	    // the curvature is too low
	    direction = 1 - direction;
	  } else {
	    if(iter >= (int)(0.5 * max_steps)) {
	      break;
	    }
	    // // check the strong Wolfe conditions
	    // if(curr_grad_curv
	    break;
	  }
	}

      } else {
	if(curr_val > best_val) {
      	  argbest_step = step;
      	  best_val = curr_val;
      	}
	if(curr_val < init_val + (control_armijo * step * init_grad_curve)) {
	  direction = 1;
	  TRACE << "insufficient progress; decreasing now";
	} else {
	  // if we've used up 3/4 of our iterations, then we'll take just satisfying the
	  // Armijo conditions
	  if(iter >= (int)(0.75 * max_steps)) {
	    break;
	  }
	  // now, what about the Wolfe (curvature) conditions?
	  double curr_grad_curv =
	    dir_correction * ferrum::dot(g_direction, grad_val);
	  if(curr_grad_curv < control_wolfe * init_grad_curve) {
	    // the curvature is too low
	    direction = 1 - direction;
	  } else {
	    if(iter >= (int)(0.5 * max_steps)) {
	      break;
	    }
	    // TODO: check the strong Wolfe conditions
	    break;
	  }
	}
      }
      
      if(step < tolerances[0]) {
	if(warn)
	  WARN << "Step = " << step << " became less than the minimum tolerance; returning step of size " << argbest_step;
	step = 0.0;
	break;
      }
      if(step > tolerances[1]) {
	if(warn) WARN << "Step = " << step << " became greater than the maximum tolerance; returning step of size " << argbest_step;
	step = 0.0;
	break;
      }
      if(iter > max_steps) {
	if(warn) WARN << "The backtracking line search is at " << iter << " steps. Something is probably wrong. Returning step of size " << argbest_step;
	step = 0.0;
	break;
      }
      //increase_step:
      step *= factors[direction];
    }
    DEBUG << "argbest_step " << argbest_step << " gives new value " << best_val << " vs. init val " << init_val;
    return argbest_step;
  }
}
