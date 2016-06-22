#ifndef CONCRETECPP_UTIL_IO_H_
#define CONCRETECPP_UTIL_IO_H_

#include <boost/log/trivial.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include "concrete/communication_types.h"

#include "thrift/protocol/TBinaryProtocol.h"
#include "thrift/protocol/TCompactProtocol.h"
#include "thrift/protocol/TProtocol.h"
#include "thrift/protocol/TProtocolException.h"

#include "thrift/transport/TTransportUtils.h"
#include "thrift/transport/TFileTransport.h"
#include "thrift/transport/TFDTransport.h"
#include "TGZipTransport.h"

#include <fcntl.h>
#include <iostream>

#include <vector>

namespace concrete { 
  namespace util {
    typedef apache::thrift::transport::TFDTransport TFDTransport;
    typedef apache::thrift::transport::TGZipTransport TZlibTransport;
    typedef apache::thrift::transport::TBufferedTransport TBufferedTransport;
    typedef apache::thrift::protocol::TBinaryProtocol TBinaryProtocol;
    typedef apache::thrift::protocol::TCompactProtocol TCompactProtocol;
    class concrete_io {
    public:
      concrete_io() {};

      template <typename P> std::vector< concrete::Communication > deserialize_compressed_communications(const char *file_name) {
      	int fd = open(file_name, O_RDONLY);

      	boost::shared_ptr<TFDTransport> innerTransport(new TFDTransport(fd));
      	boost::shared_ptr<TZlibTransport> tzt(new TZlibTransport(innerTransport));
      	boost::shared_ptr<TBufferedTransport> transport(new TBufferedTransport(tzt));
      	boost::shared_ptr<P> protocol(new P(transport));
      	transport->open();
      	std::vector< concrete::Communication > ret_list;
	bool keepReading = true;
      	while(keepReading) {
      	  concrete::Communication comm;
      	  try {
      	    comm.read(protocol.get());
	    ret_list.push_back(comm);
      	  } catch(const apache::thrift::protocol::TProtocolException &tpe) {
      	    BOOST_LOG_TRIVIAL(error) << "An exception occurred. Exception: " << tpe.what();
      	    throw tpe;
      	  } catch(const apache::thrift::transport::TTransportException &tpe) {
	    if(tpe.getType() == 3/*END_OF_FILE*/) {	      
	      keepReading = false;	      
	    } else {
	      BOOST_LOG_TRIVIAL(error) << "A TTransportException occurred. Exception: " << tpe.what();
	      throw tpe;
	    }
	  }
      	}
	try {
	  transport->close();
	} catch(const apache::thrift::transport::TTransportException& tte) {
	  BOOST_LOG_TRIVIAL(warning) << "A TTransportException has been caught in closing the GZipped compressed read-in utility. This is expected though MUST be fixed" << tte.what();
	}
      	close(fd);
      	return ret_list;
      }

      template <typename P, typename Struct> void deserialize(Struct *my_struct, const char *file_name) {
	int fd = open(file_name, O_RDONLY);

	boost::shared_ptr<TFDTransport> innerTransport(new TFDTransport(fd));
	boost::shared_ptr<TBufferedTransport> transport(new TBufferedTransport(innerTransport));
	boost::shared_ptr<P> protocol(new P(transport));
	transport->open();

	try {
	  my_struct->read(protocol.get());
	} catch(const apache::thrift::protocol::TProtocolException &tpe) {
	  BOOST_LOG_TRIVIAL(error) << "An exception occurred. Exception: " << tpe.what();
	  throw tpe;
	}
	transport->close();
	close(fd);
      };
      template <typename Struct> void deserialize_binary(Struct *my_struct, const char *file_name) {
	deserialize<TBinaryProtocol, Struct>(my_struct, file_name);
      };
      template <typename Struct> void deserialize_compact(Struct *my_struct, const char *file_name) {
	deserialize<TCompactProtocol, Struct>(my_struct, file_name);
      };
    };

    // from http://stackoverflow.com/questions/9059187/equivalent-c-to-python-generator-pattern
    class CommunicationSequence {
      typedef void (CommunicationSequence::*BoolLike)();
      void non_comparable() {
      }
    public:
      typedef concrete::Communication value_type;
    protected:
      bool done;
      value_type comm;
    public:
      typedef value_type const& reference;
      typedef value_type const* pointer;

    CommunicationSequence(): done(false) {
      }

      virtual ~CommunicationSequence() {
      }
      
      // Safe Bool idiom
      operator BoolLike() const {
    	return done ? nullptr : &CommunicationSequence::non_comparable;
      }
      reference operator*() const { return comm; }
      pointer operator->() const { return &comm; }

      virtual CommunicationSequence& begin() {
	this->operator++();
	return *this;
      }

      bool keep_reading() {
	return !done;
      }

      virtual CommunicationSequence& operator++() {
	return *this;
      }

      /* CommunicationSequence operator++(int) { */
      /* 	CommunicationSequence const tmp(*this); */
      /* 	++*this; */
      /* 	return tmp; */
      /* } */
    };

    class GZipCommunicationSequence : public CommunicationSequence {
    private:
      int fd;
      typedef TCompactProtocol P;
      boost::shared_ptr<TFDTransport> innerTransport;
      boost::shared_ptr<TZlibTransport> tzt;
      boost::shared_ptr<TBufferedTransport> transport;
      boost::shared_ptr<P> protocol;
    public:
      typedef concrete::Communication value_type;
      typedef value_type const& reference;
      typedef value_type const* pointer;

      GZipCommunicationSequence(const std::string& gzip_file_name) {
      	fd = open(gzip_file_name.c_str(), O_RDONLY);
	innerTransport = boost::shared_ptr<TFDTransport>(new TFDTransport(fd));
      	tzt = boost::shared_ptr<TZlibTransport>(new TZlibTransport(innerTransport));
      	transport = boost::shared_ptr<TBufferedTransport>(new TBufferedTransport(tzt));
      	protocol = boost::shared_ptr<P>(new P(transport));
      	transport->open();
      }

      virtual CommunicationSequence& operator++() {
    	// handle the sequence stuff here
    	assert(!done);

	concrete::Communication comm_attempt;
	try {
	  comm_attempt.read(protocol.get());
	  comm = comm_attempt;
	  return *this;
	} catch(const apache::thrift::protocol::TProtocolException &tpe) {
	  BOOST_LOG_TRIVIAL(error) << "An exception occurred. Exception: " << tpe.what();
	  throw tpe;
	} catch(const apache::thrift::transport::TTransportException &tpe) {
	  if(tpe.getType() == 3/*END_OF_FILE*/) {	      
	  } else {
	    BOOST_LOG_TRIVIAL(error) << "A TTransportException occurred. Exception: " << tpe.what();
	    throw tpe;
	  }
	}

    	//if (ij.second != Max) { ++ij.second; return *this; }
    	//if (ij.first != Max) { ij.second = 0; ++ij.first; return *this; }
    	done = true;

	try {
	  transport->close();
	} catch(const apache::thrift::transport::TTransportException& tte) {
	  BOOST_LOG_TRIVIAL(warning) << "A TTransportException has been caught in closing the GZipped compressed read-in utility. This is expected though MUST be fixed" << tte.what();
	}
      	close(fd);

    	return *this;
      }
    };

    class DirectoryCommunicationSequence : public CommunicationSequence {
    private:
      int fd;
      std::string directory_;
      boost::filesystem::path path_;
      boost::filesystem::directory_iterator directory_iterator_; 
      boost::filesystem::directory_iterator directory_end_;
      bool initialized_;
      typedef TCompactProtocol P;
      concrete_io conc_io_;

      void increment_di() {
	if(!done) {
	  ++directory_iterator_;
	  done = directory_iterator_ == directory_end_;
	}
      }

      void ensure_not_directory() {
	if(done) return;
	boost::filesystem::path f_path = boost::filesystem::path(directory_iterator_->path());
	if(boost::filesystem::is_directory(f_path) && !done) {
	  increment_di();
	  ensure_not_directory();
	}
      }

      void init() {
	initialized_ = true;
	directory_iterator_ = boost::filesystem::directory_iterator(path_);
      }
    public:
      typedef concrete::Communication value_type;
      typedef value_type const& reference;
      typedef value_type const* pointer;

    DirectoryCommunicationSequence(const std::string& directory_name) : 
      directory_(directory_name), 
	path_(directory_), initialized_(false) {
      }      

      virtual CommunicationSequence& operator++() {
	if(!initialized_) {
	  init();
	} else {
	  increment_di();
	  ensure_not_directory();
	}
    	// handle the sequence stuff here
    	if(!done) {
	  boost::filesystem::path f_path(directory_iterator_->path());
	  BOOST_LOG_TRIVIAL(debug) << "finding : " << f_path;
	  concrete::Communication comm_attempt;
	  try {
	    conc_io_.deserialize<P, concrete::Communication>(&comm_attempt, f_path.c_str());
	    comm = comm_attempt;
	    return *this;
	  } catch(const apache::thrift::protocol::TProtocolException &tpe) {
	    BOOST_LOG_TRIVIAL(error) << "An exception occurred. Exception: " << tpe.what();
	    throw tpe;
	  } catch(const apache::thrift::transport::TTransportException &tpe) {
	    if(tpe.getType() == 3/*END_OF_FILE*/) {   
	    } else {
	      BOOST_LOG_TRIVIAL(error) << "A TTransportException occurred. Exception: " << tpe.what();
	      throw tpe;
	    }
	  }
	  done = true;
	}
    	return *this;
      }
    };

    static inline void get_communication_sequence(const std::string& f_path_name, CommunicationSequence*& csp) {
      boost::filesystem::path f_path(f_path_name);
      if(boost::filesystem::is_directory(f_path)) {
	csp = new DirectoryCommunicationSequence(f_path_name);
      } else {
	csp = new GZipCommunicationSequence(f_path_name);
      }
    }
  }
}

#endif
