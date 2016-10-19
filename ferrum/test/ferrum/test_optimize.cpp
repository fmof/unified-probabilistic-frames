#ifndef PUBLIC_UPF
#include "gtest/gtest.h"

//#include "ferrum/boost/bind.hpp"
#include "ferrum/logging.hpp"
#include "ferrum/optimize.hpp"

#include <vector>

class SimpleLibLBFGSMaxent {
private:
  typedef lbfgsfloatval_t flt;
  const int counts[2] = {2, 1};

public:
  SimpleLibLBFGSMaxent() {
  }

  static flt evaluate(void *instance, const flt *x, flt *g, const int n, const flt step) {
    //int *i_counts = (int*)instance;
    int i_counts[2] = {2, 1};
    flt z = 0.;
    for (int i = 0; i < n; ++i)
      z += exp(x[i]);

    flt lz = log(z);

    flt fx = 0.;
    for (int i = 0; i < n; ++i)
      fx -= i_counts[i] * (x[i] - lz);

    for (int i = 0; i < n; ++i)
      g[i] = -i_counts[i] + 3 * exp(x[i])/z;

    return fx;
  }

  static int progress(void *instance, const flt *x, const flt *g,
		      const flt fx, const flt xnorm,
		      const flt gnorm, const flt step,
		      int n, int k, int ls) {
    return 0;
  }

  optimize::LibLBFGSFunction objective() {
    optimize::LibLBFGSFunction func;
    func.eval = &SimpleLibLBFGSMaxent::evaluate;
    func.progress = &SimpleLibLBFGSMaxent::progress;
    func.params = (void*)(this->counts);
    return func;
  }
};

TEST(liblbfgs_minimizer, maxent2_cppstyle) {
  SimpleLibLBFGSMaxent maxent;  
  optimize::LibLBFGSFunction my_func = maxent.objective();
  std::vector<double> point = {0.0, 0.0};
  optimize::LibLBFGSMinimizer optimizer(2);
  int opt_status = optimizer.minimize(&my_func, point);
  ASSERT_EQ(0, opt_status);
  ASSERT_NEAR(0.34,  point[0], 1E-2);
  ASSERT_NEAR(-0.34,  point[1], 1E-2);
}

struct Paraboloid2DParams {
  double a = 1.0;
  double x0 = 0.0;
  unsigned int e0 = 2;
  double b = 1.0;
  double y0 = 0.0;
  unsigned int e1 = 2;
  double c = 1.0;
};

/**
 * This computes
 *  y = (x-x0)^2/a^2
 */
class Parabola {
public:
  static double f(void* params, const double* point, int size) {
    Paraboloid2DParams* p = (Paraboloid2DParams*)params;
    double a = (point[0] - p->x0)/(p->a);
    double res = pow(a, p->e0);
    return res;
  }
  static void df(void* params, const double* point, double* grad, int size) {
    Paraboloid2DParams* p = (Paraboloid2DParams*)params;
    const double a2 = pow(p->a, p->e0);
    grad[0] = (p->c)*(p->e0)*pow(point[0] - p->x0, p->e0 - 1)/a2;
  }
};

TEST(AdaGrad, parabola1d_minimize_internals) {
  const size_t support = 1;
  Paraboloid2DParams p;
  p.x0 = 1;
  p.a = 2;
  optimize::AdaGrad ag(support);
  optimize::AdaGradFunction agf(&Parabola::f,
				&Parabola::df,
				(void*)&p);
  std::vector<double> point = {-10};
  ASSERT_NEAR(30.25, agf.f((void*)&p, point.data(), support), 1E-8);
  {
    std::vector<double> grad(support, 0.0);
    agf.fdf((void*)&p, point.data(), grad.data(), support);
    ASSERT_NEAR((-10.0 - 1.0)/2.0, grad[0], 1E-8);
  }
  const size_t iters = 5;
  std::vector<double> sgs(support, 0.0);
  // ferrum::print_1d(point);
  for(size_t i = 0; i < iters; ++i) {
    //std::cout << "point ";
    //ferrum::print_1d(point);
    {
      std::vector<double> grad(support, 0.0);
      agf.fdf((void*)&p, point.data(), grad.data(), support);
      //std::cout << "grad ";
      //ferrum::print_1d(grad);
      sgs[0] += grad[0]*grad[0];
    }
    ag.minimize(&agf, point, 1);
    ASSERT_NEAR(sgs[0], ag.sgs()[0], 1E-8);
  }
}

TEST(AdaGrad, parabola1d_minimize) {
  const size_t support = 1;
  Paraboloid2DParams p;
  p.x0 = 1;
  p.a = 2;
  optimize::AdaGradFunction agf(&Parabola::f,
				&Parabola::df,
				(void*)&p);
  std::vector<double> point = {-1};
  {
    optimize::AdaGrad ag(support);
    int optstat = ag.minimize(&agf, point, 1000);
    //ferrum::print_1d(point);
    const double val = agf.f((void*)&p, point.data(), support);
    INFO << optimize::Minimizer::Result::SUCCESS << " " << optimize::Minimizer::Result::MAX_ITERATIONS;
    INFO << optstat;
    if(optstat == optimize::Minimizer::Result::SUCCESS) {
      EXPECT_NEAR(0.0, val, 1E-8);
    } else {
      EXPECT_EQ(optstat, optimize::Minimizer::Result::F_EVAL_LIMIT);
    }
  }
  {
    optimize::AdaGrad ag(support);
    point[0] = -1;
    const size_t iters = 10;
    for(size_t i = 0; i < iters; ++i) {
      // std::cout << "point ";
      // ferrum::print_1d(point);
      std::vector<double> grad(support, 0.0);
      agf.fdf((void*)&p, point.data(), grad.data(), support);
      // std::cout << "grad ";
      // ferrum::print_1d(grad);
      // const double val = agf.f((void*)&p, point.data(), support);
      // std::cout << "val " << val << std::endl;
      EXPECT_EQ(optimize::Minimizer::Result::MAX_ITERATIONS,
		ag.minimize(&agf, point, 1));
    }
  }
}

/**
 * This computes
 *  z/c = (x-x0)^2/a^2 + (y-y0)^2/b^2
 */
class Paraboloid2D {
public:
  static double f(void* params, const double* point, int size) {
    Paraboloid2DParams* p = (Paraboloid2DParams*)params;
    double a = (point[0] - p->x0)/(p->a);
    double b = (point[1] - p->y0)/(p->b);
    double res = pow(a, p->e0)  + pow(b, p->e1);
    res *= (p->c);
    return res;
  }
  static void df(void* params, const double* point, double* grad, int size) {
    Paraboloid2DParams* p = (Paraboloid2DParams*)params;
    const double a2 = pow(p->a, p->e0);
    const double b2 = pow(p->b, p->e1);
    grad[0] = (p->c)*(p->e0)*pow(point[0] - p->x0, p->e0 - 1)/a2;
    grad[1] = (p->c)*(p->e1)*pow(point[1] - p->y0, p->e1 - 1)/b2;
  }
};

TEST(AdaGrad, parabola2d_minimize) {
  const size_t support = 2;
  Paraboloid2DParams p;
  p.x0 = 1;
  p.y0 = -2;
  p.a = 2;
  p.b = 4;
  optimize::AdaGrad ag(support);
  optimize::AdaGradFunction agf(&Paraboloid2D::f,
				&Paraboloid2D::df,
				(void*)&p);
  std::vector<double> point = {-10, -50};
  ASSERT_NEAR(174.25, agf.f((void*)&p, point.data(), support), 1E-8);
  {
    std::vector<double> grad(2, 0.0);
    agf.fdf((void*)&p, point.data(), grad.data(), support);
    ASSERT_NEAR((-10.0 - 1.0)/2.0, grad[0], 1E-8);
    ASSERT_NEAR((-50.0 + 2.0)/8.0, grad[1], 1E-8);
  }
  for(size_t i = 0; i < 5; ++i) {
    ag.minimize(&agf, point, 1);
    agf.f((void*)&p, point.data(), support);
    //ferrum::print_1d(point);
    std::vector<double> grad(2, 0.0);
    agf.fdf((void*)&p, point.data(), grad.data(), support);
    //ferrum::print_1d(grad);
  }
    //ag.minimize(&agf, point, 50000);
  // const double val = agf.f((void*)&p, point.data(), support);
  // ASSERT_NEAR(0.0, val, 1E-8);
  // ASSERT_NEAR(1.0, point[0], 1E-8);
  // ASSERT_NEAR(-2.0, point[1], 1E-8);
}

TEST(AdaGrad, loglin_gradop) {
  const size_t support = 2;
  const double num_obs = 21.0;
  const double obs[2] = {20.0, 1.0};
  double prob[2] = {0, 0};
  std::vector<double> point(support, 1.0);
  std::vector<double> grad(2);
  double Z = 0.0;
  for(size_t i = 0; i < support; ++i) {
    prob[i] = exp(point[i]);
    Z += prob[i];
  }
  for(size_t i = 0; i < support; ++i) {
    prob[i] /= Z;
    grad[i] = num_obs*prob[i] - obs[i];
  }
  ASSERT_NEAR(grad[0], 21.0/2.0 - 20, 1E-8);
  ASSERT_NEAR(grad[1], 21.0/2.0 - 1, 1E-8);
  optimize::AdaGrad ag(support);
  double eta = ag.eta();
  EXPECT_NEAR(0.1, eta, 1E-8);
  double eps = ag.epsilon();
  EXPECT_NEAR(1E-8, eps, 1E-10);
  std::vector<double> xpt_point = {
    point[0] - (eta / sqrt(eps + pow(grad[0],2))) * grad[0],
    point[1] - (eta / sqrt(eps + pow(grad[1],2))) * grad[1]
  };
  ag.grad_op(grad, point);
  for(size_t i = 0; i < 2; ++i) {
    ASSERT_NEAR(point[i], xpt_point[i], 1E-8);
  }
  Z = 0.0;
  for(size_t i = 0; i < support; ++i) {
    prob[i] = exp(point[i]);
    Z += prob[i];
  }
  for(size_t i = 0; i < support; ++i) {
    prob[i] /= Z;
    grad[i] = num_obs*prob[i] - obs[i];
  }
}

TEST(LibLBFGSNoOp, grad_set) {
  optimize::LibLBFGSNoOp noop;
  std::vector<double> point(4, 1.0);
  std::vector<double> grad(4, 1.0);
  noop.grad(point, grad.data());
  for(size_t k = 0; k < 4; ++k) {
    ASSERT_EQ(0.0, grad[k]);
  }
}

TEST(FunctionLinesearch, parabola2d_lower) {
  const size_t support = 2;
  Paraboloid2DParams p;
  p.x0 = 1;
  p.y0 = -2;
  p.a = 2;
  p.b = 4;
  optimize::AdaGrad ag(support);
  optimize::AdaGradFunction agf(&Paraboloid2D::f,
				&Paraboloid2D::df,
				(void*)&p);
  std::vector<double> point = {-10, -50};
  double ls = agf.linesearch(point, 1E-4, 0.9, 1.0, optimize::Function::Want::LOWER);
  INFO << "line search value is " << ls;
}

TEST(FunctionLinesearch, parabola2d_higher) {
  const size_t support = 2;
  Paraboloid2DParams p;
  p.x0 = 1;
  p.y0 = -2;
  p.a = 2;
  p.b = 4;
  p.c = -1;
  optimize::AdaGrad ag(support);
  optimize::AdaGradFunction agf(&Paraboloid2D::f,
				&Paraboloid2D::df,
				(void*)&p);
  std::vector<double> point = {-10, -50};
  double ls = agf.linesearch(point, 1E-4, 0.9, 1.0, optimize::Function::Want::HIGHER);
  INFO << "line search value is " << ls;
}
#endif
