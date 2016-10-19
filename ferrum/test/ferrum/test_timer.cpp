#include "gtest/gtest.h"

#include "ferrum/logging.hpp"
#include "ferrum/timer.hpp"


TEST(Timer, functionality) {
  {
    ferrum::Timer t(__func__);
    for(int i=0; i < 1000000; ++i) {
    }
  }
  {
    ferrum::Timer t(std::string(__func__) + ": my frame");
    for(int i=0; i < 1000000; ++i) {
    }
  }
}
