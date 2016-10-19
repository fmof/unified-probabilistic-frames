#include "gtest/gtest.h"

#include "ferrum/logging.hpp"
#include "ferrum/svi_util.hpp"

TEST(SVI, learning_rate) {
  using namespace ferrum;
  // the following defines delay_ = 1.0 and forgetting_ = 1.0
  StepSizeUpdater ssu;
  for(int i = 1; i <= 20; ++i) {
    ASSERT_NEAR(1.0/(double)i, ssu(), 1E-6);
    ++ssu;
  }
}
