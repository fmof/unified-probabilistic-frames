#include "gtest/gtest.h"

#include "ferrum/mathops.hpp"
#include "ferrum/stats.hpp"

#include <string>

TEST(Gamma, entropy) {
  ASSERT_NEAR(3.276169, ferrum::gamma::entropy(4, 3.5), 1E-5);
}

TEST(MultiNormal, sample1) {
  ferrum::ThreadRng rng;
  typedef ferrum::MultiNormal::EMat EMat;
  EMat mean;
  mean.resize(4, 1);
  mean << -10, -1, 5, 15;
  EMat covar;
  covar.setIdentity(4, 4);
  int num_samps = 100000;
  EMat ret;
  ferrum::MultiNormal::sample(num_samps, mean, covar, rng.get(), ret);
  ASSERT_EQ(ret.rows(), 4);
  ASSERT_EQ(ret.cols(), num_samps);
  EMat avgs = ret.rowwise().sum();
  avgs /= (double)num_samps;
  ASSERT_EQ(avgs.rows(), 4);
  ASSERT_EQ(avgs.cols(), 1);
  const double tol = 1E-1;
  ASSERT_NEAR(avgs(0), -10, tol);
  ASSERT_NEAR(avgs(1), -1, tol);
  ASSERT_NEAR(avgs(2), 5, tol);
  ASSERT_NEAR(avgs(3), 15, tol);
  INFO << "all good";
}
