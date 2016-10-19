#include "ferrum/logging.hpp"

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/log/attributes/current_thread_id.hpp>

#include <boost/log/core.hpp>

//This line **needs** to go here, as 
#include <boost/log/support/date_time.hpp>

#include <boost/log/expressions.hpp>

#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/named_scope.hpp>

#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/severity_logger.hpp>

#include <boost/log/trivial.hpp>

#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
  
#include <iostream>

namespace logging {
  void init() {
  }
}
#ifdef LOG_AS_COUT 
// stdout
#elif defined(LOG_AS_BOOST)
BOOST_LOG_GLOBAL_LOGGER_INIT(my_logger, logger_t) {
  logger_t lg;

  boost::log::add_common_attributes();

  boost::log::add_console_log
    (
     std::cout,
     boost::log::keywords::format =
     (
     boost::log::expressions::stream
     << "[" << boost::log::expressions::attr< unsigned int >("LineID") << "] "
     << "[" << boost::log::expressions::attr< boost::log::attributes::current_thread_id::value_type >("ThreadID") << "] "
     << "[" << boost::log::expressions::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S.%f") << "] "
     << "[" << boost::log::trivial::severity << "] "
     << boost::log::expressions::smessage
      ),
     boost::log::keywords::auto_flush = true
     );

  boost::log::core::get()->set_filter
    (
     boost::log::trivial::severity >= boost::log::trivial::info
     );

  boost::log::add_common_attributes();

  return lg;
}

#else
// GLOG
bool fLB::FLAGS_logtostderr = 0;
#endif

