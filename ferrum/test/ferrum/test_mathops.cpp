#include "gtest/gtest.h"

#include "ferrum/dmc.hpp"

#include <vector>

TEST(mathops, log_add) {
  std::vector<double> vec(3, -1.7);
  ASSERT_NEAR(-0.6013877113, mathops::log_add(vec), 1E-6);
}

TEST(mathops, sample_uniform_log) {
  double x = -.56;
  mathops::sample_uniform_log(x);
}

TEST(mathops, log_add2) {
  std::vector<double> vec(5, -1.7);
  const double norm = -0.090562087565899;
  ASSERT_NEAR(norm, mathops::log_add(vec), 1E-6);
}

TEST(mathops, log_sum_exp1) {
  std::vector<double> vec(5, -1.7);
  const double norm = -0.090562087565899;
  ASSERT_NEAR(norm, mathops::log_sum_exp(vec), 1E-6);
}

TEST(mathops, log_sum_exp2) {
  std::vector<double> vec = {-1048.45, -908.498, -1025.55, -1030.88};
  const double norm = -908.498;
  ASSERT_EQ(-908.498, ferrum::max(vec));
  std::vector<double> exped;
  std::vector<double> expected_exped = {1.6581303503252086e-61, 1.0, 1.4620502664253475e-51, 7.082273851990543e-54};
  for(size_t i = 0; i < 4; ++i) {
    double e = mathops::exp( vec[i] - ferrum::max(vec));
    exped.push_back(e);
    EXPECT_NEAR(expected_exped[i], e, 1E-3) << " coordinate " << i << " failed";
  }
  ASSERT_NEAR(norm, mathops::log_sum_exp(vec), 1E-3);
}

TEST(ferrum, prob_from_unnorm_lp) {
  std::vector<double> vec = {-1048.45, -908.498, -1025.55, -1030.88};
  const double norm = -908.498;
  ASSERT_NEAR(norm, mathops::log_sum_exp(vec), 1E-3);
  std::vector<double> vec2(vec);
  ferrum::sum(-1 * norm, &vec2);
  std::vector<double> foo = {-139.952, 0.0, -117.05199999999991, -122.38200000000006};
  for(size_t i = 0; i < 4; ++i) {
    EXPECT_NEAR(foo[i], vec2[i], 1E-3) << " coordinate " << i << " failed";
  }
  ferrum::exp(&vec2);
  ferrum::prob_from_unnorm_lp(&vec);
  std::vector<double> expected_exped = {1.6581303503252086e-61, 1.0, 1.4620502664253475e-51, 7.082273851990543e-54};
  for(size_t i = 0; i < 4; ++i) {
    EXPECT_NEAR(expected_exped[i], vec2[i], 1E-3) << " coordinate " << i << " failed";
    EXPECT_NEAR(expected_exped[i], vec[i], 1E-3) << " coordinate " << i << " failed";
  }
}
