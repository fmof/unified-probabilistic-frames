#ifndef FERRUM_LIBNAR_TIMER_HPP_
#define FERRUM_LIBNAR_TIMER_HPP_

#include <chrono>
#include <ctime>
#include <string>

namespace ferrum {
  class Timer {
  public:
    Timer(const std::string& frame);
    ~Timer();
    void snapshot();
    typedef std::chrono::duration<double> Duration;
    Duration lap();
  private:
    std::string frame_;
    std::chrono::time_point<std::chrono::system_clock> start_;
    std::chrono::time_point<std::chrono::system_clock> end_;
    std::chrono::time_point<std::chrono::system_clock> lap_start_;
  };
}

#endif
