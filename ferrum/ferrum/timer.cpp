#include "ferrum/logging.hpp"
#include "ferrum/timer.hpp"

namespace ferrum {
  Timer::Timer(const std::string& frame) :
    frame_(frame),
    start_(std::chrono::system_clock::now()),
    lap_start_(start_) {
  }
  Timer::~Timer() {
    end_ = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_secs = end_ - start_;
    INFO << "Elapsed time of " << elapsed_secs.count() << "s for frame " << frame_;
  }
  void Timer::snapshot() {
    std::chrono::time_point<std::chrono::system_clock> pend = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_secs = pend - start_;
    INFO << "(partial) Elapsed time of " << elapsed_secs.count() << "s for frame " << frame_;
  }
  typename Timer::Duration Timer::lap() {
    std::chrono::time_point<std::chrono::system_clock> pend = std::chrono::system_clock::now();
    Duration elap = pend - lap_start_;
    lap_start_ = pend;
    return elap;
  }
}
