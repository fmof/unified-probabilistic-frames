#ifndef PUBLIC_UPF
#include "gtest/gtest.h"

#include <boost/tokenizer.hpp>
#include "ferrum/cblas_cpp.hpp"
#include "ferrum/loglin.hpp"
#include "ferrum/mathops.hpp"
#include "ferrum/optimize.hpp"
#include "ferrum/util.hpp"
//#include "gsl/mathops::exp.h"
#include "gsl/gsl_sf_log.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace loglin;

typedef UnigramMaxent<size_t, std::vector<double> > Model;

TEST(loglin, create_model) {  
  Model ment;
  //std::vector<double> vec(3, -1.7);
  //ASSERT_NEAR(-0.6013877113, mathops::log_add(vec), 1E-6);
}

TEST(loglin, copy_assign) {
  Model m1;
  std::vector<double> w(2);
  w[0] = .5;
  w[1] = -.1;
  m1.weights(w);
  Model m2 = m1;
}

TEST(loglin, copy_ctor) {
  Model m1;
  std::vector<double> w(2);
  w[0] = .5;
  w[1] = -.1;
  m1.weights(w);
  Model m2(m1);
}

TEST(loglin, move_assign) {
  Model m1;
  std::vector<double> w(2);
  w[0] = .5;
  w[1] = -.1;
  m1.weights(w);
  Model m2 = std::move(m1);
  ASSERT_TRUE(m1.weights() == NULL);
}

TEST(loglin, move_ctor) {
  Model m1;
  std::vector<double> w(2);
  w[0] = .5;
  w[1] = -.5;
  m1.weights(w);
  Model m2(std::move(m1));
  ASSERT_TRUE(m1.weights() == NULL);
  ASSERT_NEAR(0.731059, m2.p(0), 1E-6);
}

TEST(loglin, lesson1) {
  loglin::UnigramMaxent<size_t, std::vector<double> > model;
  std::vector<double> weights(4);
  weights[0] = 1; // striped circle
  weights[1] = 1; // solid triangle
  weights[2] = 2; // solid circle
  weights[3] = 0; // striped triangle
  const double test_tol = 1E-6;
  // compute normalization
  model.weights(weights);
  const double log_norm = model.log_normalizer();
  const double expected_log_norm = 2.626523375;
  ASSERT_NEAR(expected_log_norm, log_norm, test_tol);
  // test individual log probs
  ASSERT_NEAR(1 - expected_log_norm, model.lp(0), test_tol);
  ASSERT_NEAR(1 - expected_log_norm, model.lp(1), test_tol);
  ASSERT_NEAR(2 - expected_log_norm, model.lp(2), test_tol);
  ASSERT_NEAR(0 - expected_log_norm, model.lp(3), test_tol);
  // test log likelihood
  std::map<size_t, int> data;
  data[0] = 15;
  data[1] = 10;
  data[2] = 30;
  data[3] = 5;
  const double expected_ll = -72.59140250218672;
  ASSERT_NEAR(expected_ll, model.ll(data), test_tol);
  // test gradient of ll
  const std::vector<double> grad = model.ll_grad(data);
  EXPECT_NEAR(30 - 60*model.p(2), grad[2], test_tol);
  EXPECT_NEAR(10 - 60*model.p(1), grad[1], test_tol);
  EXPECT_NEAR(15 - 60*model.p(0), grad[0], test_tol);
  EXPECT_NEAR(5  - 60*model.p(3), grad[3], test_tol);
}

struct UnigramClosure {
  typedef loglin::UnigramMaxent<int, std::vector<double> > MaxentModel;
  MaxentModel* model;
  std::vector<double>* counts;
};

class TestUnigramMaxentLM {
public:
  typedef UnigramClosure ClosureType;

  static double
  lbfgs_elbo(void *fparams, const lbfgsfloatval_t *trial_weights,
	     lbfgsfloatval_t *lbfgs_grad, const int n,
	     const lbfgsfloatval_t step) {
    ClosureType* closure = (ClosureType*)fparams;
    std::vector<double> nweights = 
      optimize::LibLBFGSVector::to_container<std::vector<double> >(trial_weights, n);
    // reset the maxent weights; this will renormalize everything
    closure->model->weights(nweights);
    double ll = closure->model->ll_dense_data(*(closure->counts));
    typedef std::vector<double> WeightType;
    WeightType grad = closure->model->ll_grad_dense_data(*(closure->counts));
    // negate
    ferrum::scalar_product(-1.0, &grad);
    optimize::LibLBFGSVector::copy(grad, lbfgs_grad);
    return -ll;
  }

  // static double
  // elbo(void *fparams, const double *trial_weights,
  //      double *grad, int n) {
  //   ClosureType* closure = (ClosureType*)fparams;
  //   std::vector<double> nweights = 
  //     optimize::LibLBFGSVector::to_container<std::vector<double> >(trial_weights, n);
  //   // reset the maxent weights; this will renormalize everything
  //   closure->model->weights(nweights);

  // }
  static double
  elbo_eval(void *fparams, const double *trial_weights, int n) {
    ClosureType* closure = (ClosureType*)fparams;
    std::vector<double> nweights = 
      optimize::LibLBFGSVector::to_container<std::vector<double> >(trial_weights, n);
    // reset the maxent weights; this will renormalize everything
    closure->model->weights(nweights);
    double ll = closure->model->ll_dense_data(*(closure->counts));
    return -ll;
  }
  static void
  elbo_grad(void *fparams, const double *trial_weights,
	    double *grad, int n) {
    ClosureType* closure = (ClosureType*)fparams;
    std::vector<double> nweights = 
      optimize::LibLBFGSVector::to_container<std::vector<double> >(trial_weights, n);
    // reset the maxent weights; this will renormalize everything
    closure->model->weights(nweights);
    closure->model->ll_grad_dense_data(*(closure->counts), grad);
    // negate
    cblas_dscal(n, -1.0, grad, 1);
  }
  optimize::LibLBFGSFunction get_liblbfgs_fdf(ClosureType* params, const size_t size) {
    optimize::LibLBFGSFunction objective;
    objective.eval   = &TestUnigramMaxentLM::lbfgs_elbo;
    objective.progress = &optimize::LibLBFGSNoOp::progress;
    objective.params = (void*)params;
    return objective;
  }

  std::shared_ptr<optimize::Function>
  get_optimizable_func(ClosureType* params, optimize::OptimizationMethod opt_how) {
    std::shared_ptr<optimize::Function> func;
    if(opt_how == optimize::OptimizationMethod::LBFGS) {
      func = std::shared_ptr<optimize::Function>
	(new optimize::LibLBFGSFunction(&TestUnigramMaxentLM::lbfgs_elbo,
					&optimize::LibLBFGSNoOp::progress,
					(void*)params));
    }
    else {
      switch(opt_how) {
      case optimize::OptimizationMethod::ADAGRAD:
	{
	  optimize::MinimizationFunctionCreator<optimize::OptimizationMethod::ADAGRAD> mfc;
	  func = mfc(&TestUnigramMaxentLM::elbo_eval,
		     &TestUnigramMaxentLM::elbo_grad,
		     (void*)params);
	}
	break;
      default:
	WARN << "Trying to create function for optimization method " << opt_how << "; unexpected, but it might work";
	break;
      }
    }
    return func;
  }

  std::vector<double> get_optimization_initial_point(int size_) {
    std::vector<double> foo(size_, 0.0);
    return foo;
  }
};

TEST(loglin, word_unigram_lm_2words_liblbfgs) {
  loglin::UnigramMaxent<int, std::vector<double> > model;
  std::map<std::string, int> vocab;
  std::vector<double> counts;
  std::string text = "the the man";
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
  int tok_count = 0;
  for (const auto& t : tokens) {
    std::map<std::string,int>::iterator it = vocab.find(t);
    if(it == vocab.end()) {
      vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
      counts.push_back(1.0);
    } else {
      counts[it->second]++;
    }    
  }
  // Test basic data stats 
  {
    ASSERT_EQ(2, vocab.size());
    ASSERT_EQ(2, counts[0]); // the
    ASSERT_EQ(1, counts[1]); // man
  }
  UnigramClosure params;
  params.counts = &counts;
  params.model = &model;
  TestUnigramMaxentLM lm;
  optimize::LibLBFGSFunction my_func = lm.get_liblbfgs_fdf(&params, vocab.size());
  std::vector<double> point = lm.get_optimization_initial_point(vocab.size());
  optimize::LibLBFGSMinimizer optimizer(point.size());
  int opt_status = optimizer.minimize(&my_func, point);
  ASSERT_EQ(0, opt_status);
  ASSERT_NEAR(point[0] - point[1], gsl_sf_log(2), 1E-6);
}

TEST(loglin, word_unigram_lm_2words_ll_grad_dense_data) {
  loglin::UnigramMaxent<int, std::vector<double> > model;
  std::map<std::string, int> vocab;
  std::vector<double> counts;
  std::string text = "the the man";
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
  int tok_count = 0;
  for (const auto& t : tokens) {
    std::map<std::string,int>::iterator it = vocab.find(t);
    if(it == vocab.end()) {
      vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
      counts.push_back(1.0);
    } else {
      counts[it->second]++;
    }    
  }
  // Test basic data stats 
  {
    ASSERT_EQ(2, vocab.size());
    ASSERT_EQ(2, counts[0]); // the
    ASSERT_EQ(1, counts[1]); // man
  }
  std::vector<double> point(vocab.size(), 0.0);
  UnigramClosure params;
  params.counts = &counts;
  params.model = &model;
  std::vector<double> grad(vocab.size(), 0.0);
  TestUnigramMaxentLM::elbo_grad((void*)(&params),
				 point.data(),
				 grad.data(),
				 vocab.size());
  std::vector<double> expected_obs = {2, 1};
  std::vector<double> expected_expt = 
    {
      3.0 * 1.0/2.0,
      3.0 * 1.0/2.0
    };
  for(size_t i = 0; i < vocab.size(); ++i) {
    EXPECT_NEAR(grad[i], expected_expt[i] - expected_obs[i], 1E-6);
  }
}

TEST(loglin, word_unigram_lm_2words_liblbfgs_wrapper) {
  loglin::UnigramMaxent<int, std::vector<double> > model;
  std::map<std::string, int> vocab;
  std::vector<double> counts;
  std::string text = "the the man";
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
  int tok_count = 0;
  for (const auto& t : tokens) {
    std::map<std::string,int>::iterator it = vocab.find(t);
    if(it == vocab.end()) {
      vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
      counts.push_back(1.0);
    } else {
      counts[it->second]++;
    }    
  }
  // Test basic data stats 
  {
    ASSERT_EQ(2, vocab.size());
    ASSERT_EQ(2, counts[0]); // the
    ASSERT_EQ(1, counts[1]); // man
  }
  UnigramClosure params;
  params.counts = &counts;
  params.model = &model;
  TestUnigramMaxentLM lm;
  optimize::OptimizationMethod how = optimize::OptimizationMethod::LBFGS;
  std::shared_ptr<optimize::Function> func = 
    lm.get_optimizable_func(&params, how);
  std::vector<double> point = lm.get_optimization_initial_point(vocab.size());
  //optimize::LibLBFGSMinimizer optimizer(point.size());
  std::shared_ptr<optimize::Minimizer> minimizer =
    optimize::Minimizer::make(how, vocab.size());
  int opt_status = minimizer->minimize(func.get(), point);
  ASSERT_EQ(0, opt_status);
  ASSERT_NEAR(point[0] - point[1], gsl_sf_log(2), 1E-6);
}

TEST(loglin, word_unigram_lm_2words_adagrad_wrapper) {
  loglin::UnigramMaxent<int, std::vector<double> > model;
  std::map<std::string, int> vocab;
  std::vector<double> counts;
  std::string text = "the the man";
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
  int tok_count = 0;
  for (const auto& t : tokens) {
    std::map<std::string,int>::iterator it = vocab.find(t);
    if(it == vocab.end()) {
      vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
      counts.push_back(1.0);
    } else {
      counts[it->second]++;
    }    
  }
  // Test basic data stats 
  {
    ASSERT_EQ(2, vocab.size());
    ASSERT_EQ(2, counts[0]); // the
    ASSERT_EQ(1, counts[1]); // man
  }
  UnigramClosure params;
  params.counts = &counts;
  params.model = &model;
  TestUnigramMaxentLM lm;
  std::vector<double> point = lm.get_optimization_initial_point(vocab.size());
  double init_val = TestUnigramMaxentLM::elbo_eval((void*)&params,
						   point.data(),
						   vocab.size());
  optimize::OptimizationMethod how = optimize::OptimizationMethod::ADAGRAD;
  std::shared_ptr<optimize::Function> func = 
    lm.get_optimizable_func(&params, how);
  std::shared_ptr<optimize::Minimizer> minimizer =
    optimize::Minimizer::make(how, vocab.size());
  int opt_status = minimizer->minimize(func.get(), point);
  ASSERT_EQ(optimize::Minimizer::Result::MAX_ITERATIONS, opt_status);
  opt_status = dynamic_cast<optimize::AdaGrad*>(minimizer.get())->minimize(func.get(), point, 300);
  ASSERT_EQ(optimize::Minimizer::Result::F_EVAL_LIMIT, opt_status);
  double final_val = TestUnigramMaxentLM::elbo_eval((void*)&params,
						    point.data(),
						    vocab.size());
  EXPECT_LE(final_val, init_val) << "The initial value " << init_val << " is not greater than the final value " << final_val;
  //ASSERT_NEAR(point[0] - point[1], gsl_sf_log(2), 1E-6) << "Failed on " << how;
}

// TEST(loglin, word_unigram_lm_2words_liblbfgs_adagrad_wrapper) {
//   loglin::UnigramMaxent<int, std::vector<double> > model;
//   std::map<std::string, int> vocab;
//   std::vector<double> counts;
//   std::string text = "the the man";
//   boost::char_separator<char> sep(" ");
//   boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
//   int tok_count = 0;
//   for (const auto& t : tokens) {
//     std::map<std::string,int>::iterator it = vocab.find(t);
//     if(it == vocab.end()) {
//       vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
//       counts.push_back(1.0);
//     } else {
//       counts[it->second]++;
//     }    
//   }
//   // Test basic data stats 
//   {
//     ASSERT_EQ(2, vocab.size());
//     ASSERT_EQ(2, counts[0]); // the
//     ASSERT_EQ(1, counts[1]); // man
//   }
//   UnigramClosure params;
//   params.counts = &counts;
//   params.model = &model;
//   // for(optimize::OptimizationMethod how : {optimize::OptimizationMethod::LBFGS,
//   // 	optimize::OptimizationMethod::ADAGRAD} ) {
//   TestUnigramMaxentLM lm;
//   std::vector<double> point = lm.get_optimization_initial_point(vocab.size());
//   // optimize::OptimizationMethod how = optimize::OptimizationMethod::LBFGS;
//   // std::shared_ptr<optimize::Function> func = 
//   //   lm.get_optimizable_func(&params, how);
//   // std::shared_ptr<optimize::Minimizer> minimizer =
//   //   optimize::Minimizer::make(how, vocab.size());
//   // int opt_status = minimizer->minimize(func.get(), point);
//   // ASSERT_EQ(0, opt_status);
//   // ASSERT_NEAR(point[0] - point[1], gsl_sf_log(2), 1E-6) << "Failed on " << how;
//     //}
// }


TEST(loglin, ll_grad_dense_data) {
  loglin::UnigramMaxent<int, std::vector<double> > model;
  std::map<std::string, int> vocab;
  std::vector<double> counts;
  std::string text = "the the man";
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(text, sep);
  int tok_count = 0;
  for (const auto& t : tokens) {
    std::map<std::string,int>::iterator it = vocab.find(t);
    if(it == vocab.end()) {
      vocab.insert(it, std::pair<std::string, int>(t, tok_count++));
      counts.push_back(1.0);
    } else {
      counts[it->second]++;
    }    
  }
  // Test basic data stats 
  {
    ASSERT_EQ(2, vocab.size());
    ASSERT_EQ(2, counts[0]); // the
    ASSERT_EQ(1, counts[1]); // man
  }
  UnigramClosure params;
  params.counts = &counts;
  params.model = &model;
  TestUnigramMaxentLM lm;
  std::vector<double> point = lm.get_optimization_initial_point(vocab.size());
  model.weights(point);
  typedef std::vector<double> WeightType;
  WeightType grad1 = model.ll_grad_dense_data(*(params.counts));
  WeightType grad2(grad1.size(), 0.0);
  model.ll_grad_dense_data(*params.counts, grad2.data());
  for(size_t i = 0; i < grad1.size(); ++i) {
    ASSERT_NEAR(grad1[i], grad2[i], 1E-8);
  }
}
#endif
