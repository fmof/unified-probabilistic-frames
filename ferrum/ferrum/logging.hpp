#ifndef FERRUM_LOGGING_H_
#define FERRUM_LOGGING_H_

namespace logging {  
  void init();
}

#ifdef LOG_AS_COUT
#include <iostream>

#define ADD_FILE_LINE_FUNC << "(" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << ") "

#define TRACE std::cout << std::endl << "TRACE " << ADD_FILE_LINE_FUNC
#define DEBUG std::cout << std::endl << "DEBUG " << ADD_FILE_LINE_FUNC
#define INFO std::cout << std::endl << "INFO " << ADD_FILE_LINE_FUNC
#define WARN std::cout << std::endl << "WARN  " << ADD_FILE_LINE_FUNC
#define ERROR std::cout << std::endl << "ERROR "  << ADD_FILE_LINE_FUNC

#elif defined(LOG_AS_BOOST)

#include <boost/log/core.hpp>

#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/severity_logger.hpp>


#define ADD_FILE_LINE_FUNC << "[" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "] "

#define TRACE BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::trace) 
#define DEBUG BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::debug) 
#define INFO BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::info) 
#define WARN BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::warning) 
#define ERROR BOOST_LOG_SEV(my_logger::get(), boost::log::trivial::error) 

//Narrow-char thread-safe logger.
typedef boost::log::sources::severity_logger_mt< boost::log::trivial::severity_level > logger_t;

//declares a global logger with a custom initialization
// behind the scenes, this macro defines a struct
BOOST_LOG_GLOBAL_LOGGER(my_logger, logger_t)

// #define TRACE_TRIVIAL TRACE ADD_FILE_LINE_FUNC
// #define DEBUG_TRIVIAL DEBUG ADD_FILE_LINE_FUNC
// #define INFO_TRIVIAL INFO ADD_FILE_LINE_FUNC
// #define WARN_TRIVIAL WARN ADD_FILE_LINE_FUNC
// #define ERROR_TRIVIAL ERROR ADD_FILE_LINE_FUNC

// #define TRACE_CUSTOM TRACE 
// #define DEBUG_CUSTOM DEBUG 
// #define INFO_CUSTOM INFO 
// #define WARN_CUSTOM WARN 
// #define ERROR_CUSTOM ERROR 

// #define GET_MACRO(_0, _1, _2, NAME, ...) NAME

// #define TRACE(...) GET_MACRO(_0, ##__VA_ARGS__, TRACE_TRIVIAL, TRACE_CUSTOM)(__VA_ARGS__)
// #define DEBUG(...) GET_MACRO(_0, ##__VA_ARGS__, TRACE_TRIVIAL, TRACE_CUSTOM)(__VA_ARGS__)
// #define INFO(...) GET_MACRO(_0, ##__VA_ARGS__, TRACE_TRIVIAL, TRACE_CUSTOM)(__VA_ARGS__)
// #define WARN(...) GET_MACRO(_0, ##__VA_ARGS__, TRACE_TRIVIAL, TRACE_CUSTOM)(__VA_ARGS__)
// #define ERROR(...) GET_MACRO(_0, ##__VA_ARGS__, TRACE_TRIVIAL, TRACE_CUSTOM)(__VA_ARGS__)

#else

#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>

#define ADD_FILE_LINE_FUNC << "[" << __FILE__ << ":" << __LINE__ << ":" << __FUNCTION__ << "] "

#if defined(FERRUM_DEBUG)
#define DEBUG COMPACT_GOOGLE_LOG_INFO.stream()
#if defined(FERRUM_UNIT_TEST) || defined(IET_UNIT_TEST)
#define TRACE COMPACT_GOOGLE_LOG_INFO.stream()
#else
#define TRACE google::NullStream().stream()
#endif
#elif defined(NDEBUG) || !defined(DEBUG)
#define TRACE google::NullStream().stream()
#define DEBUG google::NullStream().stream()
#else
#define TRACE COMPACT_GOOGLE_LOG_INFO.stream()
#define DEBUG COMPACT_GOOGLE_LOG_INFO.stream()
#endif

#define INFO COMPACT_GOOGLE_LOG_INFO.stream()
#define WARN COMPACT_GOOGLE_LOG_WARNING.stream()
#define ERROR COMPACT_GOOGLE_LOG_ERROR.stream()

#endif

#endif

