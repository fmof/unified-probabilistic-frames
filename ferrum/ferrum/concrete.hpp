/**
 * A libnar-specific utility library for Concrete.
 */

#ifndef FERRUM_LIBNAR_CONCRETE_H_
#define FERRUM_LIBNAR_CONCRETE_H_

#include "concrete/communication_types.h"
#include "concrete/entities_types.h"
#include "concrete/metadata_types.h"
#include "concrete/situations_types.h"
#include "concrete/spans_types.h"
#include "concrete/structure_types.h"
#include "concrete/uuid_types.h"
#include "concrete/version.hpp"

#include "ferrum/data_util.hpp"
#include "ferrum/util.hpp"
#include "ferrum/thrift_protocol_defs.hpp"
#include "ferrum/thrift_smart_writer.hpp"

#include <array>
#include <list>
#include <map>
#include <tuple>
#include <vector>
#include <unordered_map>

#include "concrete_util/uuid_util.h"
#include "concrete_util/io.h"
#include <chrono>
#include <fcntl.h>
#include <iostream>
#include <thrift/protocol/TJSONProtocol.h>
#include <vector>

namespace concrete {
  namespace util {
    typedef apache::thrift::transport::TFDTransport TFDTransport;
    typedef apache::thrift::transport::TGZipTransport TZlibTransport;
    typedef apache::thrift::transport::TBufferedTransport TBufferedTransport;
    typedef ferrum::thrift::TBinaryProtocol TBinaryProtocol;
    typedef ferrum::thrift::TCompactProtocol TCompactProtocol;
    typedef ferrum::thrift::TJSONProtocol TJSONProtocol;

    class TBinaryGZipCommunicationSequence : public CommunicationSequence {
    private:
      int fd;
      typedef TBinaryProtocol P;
      boost::shared_ptr<TFDTransport> innerTransport;
      boost::shared_ptr<TZlibTransport> tzt;
      boost::shared_ptr<TBufferedTransport> transport;
      boost::shared_ptr<P> protocol;
    public:
      typedef concrete::Communication value_type;
      typedef value_type const& reference;
      typedef value_type const* pointer;

      TBinaryGZipCommunicationSequence(const std::string& gzip_file_name) {
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
	  ERROR << "An exception occurred. Exception: " << tpe.what();
	  throw tpe;
	} catch(const apache::thrift::transport::TTransportException &tpe) {
	  if(tpe.getType() == 3/*END_OF_FILE*/) {	      
	  } else {
	    ERROR << "A TTransportException occurred. Exception: " << tpe.what();
	    throw tpe;
	  }
	}

    	//if (ij.second != Max) { ++ij.second; return *this; }
    	//if (ij.first != Max) { ij.second = 0; ++ij.first; return *this; }
    	done = true;

	try {
	  transport->close();
	} catch(const apache::thrift::transport::TTransportException& tte) {
	  WARN << "A TTransportException has been caught in closing the GZipped compressed read-in utility. This is expected though MUST be fixed" << tte.what();
	}
      	close(fd);

    	return *this;
      }
    };

    class TBinaryDirectoryCommunicationSequence : public CommunicationSequence {
    private:
      std::string directory_;
      boost::filesystem::path path_;
      boost::filesystem::directory_iterator directory_iterator_; 
      boost::filesystem::directory_iterator directory_end_;
      bool initialized_;
      typedef TBinaryProtocol P;
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

    TBinaryDirectoryCommunicationSequence(const std::string& directory_name) : 
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
	  DEBUG << "finding : " << f_path;
	  concrete::Communication comm_attempt;
	  try {
	    conc_io_.deserialize<P, concrete::Communication>(&comm_attempt, f_path.c_str());
	    comm = comm_attempt;
	    return *this;
	  } catch(const apache::thrift::protocol::TProtocolException &tpe) {
	    ERROR << "An exception occurred. Exception: " << tpe.what();
	    throw tpe;
	  } catch(const apache::thrift::transport::TTransportException &tpe) {
	    if(tpe.getType() == 3/*END_OF_FILE*/) {   
	    } else {
	      ERROR << "A TTransportException occurred. Exception: " << tpe.what();
	      throw tpe;
	    }
	  }
	  done = true;
	}
    	return *this;
      }
    };

    void get_binary_communication_sequence(const std::string& f_path_name, CommunicationSequence*& csp);
    //void get_compact_communication_sequence(const std::string& f_path_name, CommunicationSequence*& csp);

    /**
     * return true iff first is fully-contained in second
     */
    bool contained_in(const concrete::TextSpan& first, const concrete::TextSpan& second);

    enum ConcreteReadResult {
      OK = 0,
      E_BAD_SIZE = 1,
      E_OOM = 2,
    };

    struct ConcreteResultStruct {
      ConcreteReadResult crr;
      size_t size;
    };

    template <typename Archiver, typename ConcreteObj, typename Loader, typename... ArchiverSizeArgs>
    ConcreteResultStruct load_from_archive(Archiver& archive, ConcreteObj& obj, Loader loader, ArchiverSizeArgs& ... size_args) {
      ConcreteResultStruct crs;
      crs.size = archive.size(size_args...);
      void* buffer = malloc(crs.size);
      if(! buffer) {
	crs.crr = ConcreteReadResult::E_OOM;
      } else {
	TRACE << "sizes: " << crs.size;
	archive.read(buffer, crs.size);
	if(crs.size == 0) {
	  crs.crr = ConcreteReadResult::E_BAD_SIZE;
	} {
	  loader//(
		 //ferrum::thrift::thrift_struct_from_buffer<Protocol>
	    (
	     buffer,
	     crs.size,
	     &obj
	     );
	  free(buffer);
	  crs.crr = ConcreteReadResult::OK;
	}
      }
      return crs;
    }

    template <typename P>
    class ConcreteSmartWriter : public ferrum::thrift::ThriftSmartWriter<P> {
    private:
      int fd_;
      boost::shared_ptr<TFDTransport> innerTransport_;
      boost::shared_ptr<TBufferedTransport> transport_;
    public:
      ConcreteSmartWriter<P>(const std::string& base_name);
      virtual ~ConcreteSmartWriter();
      P* get();
      P* get(const std::string& suffix);
      int close_csw();
      virtual ConcreteSmartWriter<P>* clone() const;
    protected:
      void _load(void* sv);
      void _save();
    };

    template <typename CMention>
    inline bool index_overlap_mention(const CMention& mention, int idx) {
      bool found = false;
      int anchor = -1;
      if((anchor = mention.tokens.anchorTokenIndex) >= 0 &&
	 anchor == idx) {
	found = true;
      } else { // default to looking through the token index list	    
	for(int tok_idx : mention.tokens.tokenIndexList) {
	  found = (tok_idx == idx);
	  if(found) break;
	}
      }
      return found;
    }

    concrete::AnnotationMetadata make_metadata(const std::string& tool);
    template <typename C>
    inline void add_uuid(C& obj, concrete::util::uuid_factory& uf) {
      concrete::UUID uuid;
      uuid.__set_uuidString(uf.get_uuid());
      obj.__set_uuid(uuid);
    }
    template <typename C>
    inline void add_metadata(C& obj, const std::string& tool) {
      concrete::AnnotationMetadata am =
	concrete::util::make_metadata(tool);
      obj.__set_metadata(am);
    }

    // // two forward-declares
    // template <typename C, typename... ITs>   
    // bool operator==(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1,
    // 		    const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i2);
    // template <typename C, typename... ITs>
    // bool operator!=(const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i1,
    // 		    const typename StructIterator<C,ITs...>::ConcreteStructIteratorReturn& i2);

    template <typename C, typename... IterTypes>
    class StructIterator {
    protected:

      // template <int... Is>
      // struct index {
      // };
      // template <int N, int... Is>
      // struct sequence : sequence<N - 1, N - 1, Is...> {
      // };
      // template <int... Is>
      // struct sequence<0, Is...> : index<Is...> {
      // };
      template <size_t> struct index { };
      template <typename Iter>
      struct IterStruct {
	IterStruct();
    	Iter c; // current value
    	Iter e; // end point
      };
      
      class ConcreteStructIteratorReturn {
      public:
    	//ConcreteStructIteratorReturn(const concrete::Communication&);
	ConcreteStructIteratorReturn(StructIterator<C,IterTypes...>*);
    	typedef const C& reference;
    	typedef const C* pointer;
    	pointer operator->();
    	reference operator*();
	/**
	 * Advance the iterator forward; guaranteed to point to a
	 * (de)referencable address, or `ConcreteStructIteratorReturn::end()`.
	 */
    	ConcreteStructIteratorReturn& operator++(); // prefix
    	static ConcreteStructIteratorReturn end();
    	bool operator==(const typename StructIterator<C,IterTypes...>::ConcreteStructIteratorReturn&) const;
	bool operator!=(const typename StructIterator<C,IterTypes...>::ConcreteStructIteratorReturn&) const;
      private:
    	const C* obj_ptr_;
	StructIterator<C,IterTypes...>* par_;
	std::tuple< IterStruct<IterTypes>... > iters_;
	std::array< bool, sizeof...(IterTypes) > need_to_reset_;
	/**
	 * Initialize the iterators, traversing outer-to-inner indices.
	 */
	void init_iters_();
    	ConcreteStructIteratorReturn& advance_(bool, bool);
	template <size_t i> void reset_inner_(index<i>);
	void set_obj_ptr_();
	template <int> int check(bool);
      protected:
    	ConcreteStructIteratorReturn();
      };

    public:
      StructIterator(const concrete::Communication&);
      typedef const C& reference;
      typedef C value_type;
      typedef const C* pointer;
      typedef ConcreteStructIteratorReturn iterator;
      typedef const ConcreteStructIteratorReturn const_iterator;
      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;
    private:
      const concrete::Communication* c_;
    };

    template <>
    class StructIterator<concrete::Communication> : public StructIterator<concrete::Communication, int> {
    public:
      using StructIterator<concrete::Communication, int>::StructIterator;
    };
    template <>
    class StructIterator<concrete::Tokenization>  : public StructIterator<concrete::Tokenization,	std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator> {
    public:
      using StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator>::StructIterator;
    };
    template <>
    class StructIterator<Token> : public StructIterator<concrete::Token,
							std::vector<concrete::Section>::const_iterator,
							std::vector<concrete::Sentence>::const_iterator,
							std::vector<concrete::Token>::const_iterator> {
    public:
      using StructIterator<concrete::Token,
			   std::vector<concrete::Section>::const_iterator,
			   std::vector<concrete::Sentence>::const_iterator,
			   std::vector<concrete::Token>::const_iterator>::StructIterator;
    };

    // template <>
    // class StructIterator<TaggedToken> : public StructIterator<concrete::TaggedToken,
    // 							std::vector<concrete::Section>::const_iterator,
    // 							std::vector<concrete::Sentence>::const_iterator,
    // 							std::vector<concrete::TaggedToken>::const_iterator> {
    // public:
    //   using StructIterator<concrete::TaggedToken,
    // 			   std::vector<concrete::Section>::const_iterator,
    // 			   std::vector<concrete::Sentence>::const_iterator,
    // 			   std::vector<concrete::TaggedToken>::const_iterator>::StructIterator;
    // };

    class CommunicationIterator : public StructIterator<concrete::Communication, int> {
    public:
      using StructIterator<concrete::Communication, int>::StructIterator;
    };
    class TokenizationIterator  : public StructIterator<concrete::Tokenization,	std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator> {
    public:
      using StructIterator<concrete::Tokenization, std::vector<concrete::Section>::const_iterator, std::vector<concrete::Sentence>::const_iterator>::StructIterator;
    };
    class TokenIterator : public StructIterator<concrete::Token,
						std::vector<concrete::Section>::const_iterator,
						std::vector<concrete::Sentence>::const_iterator,
						std::vector<concrete::Token>::const_iterator> {
    public:
      using StructIterator<concrete::Token,
						std::vector<concrete::Section>::const_iterator,
						std::vector<concrete::Sentence>::const_iterator,
			   std::vector<concrete::Token>::const_iterator>::StructIterator;
    };


  } // end namespace util
} // end namespace concrete

namespace concrete_util {
  struct uuid_hash { 
    const size_t operator()(const concrete::UUID& uuid) const;
  };

  template <typename T > using uuid_map = std::unordered_map< const concrete::UUID, T , concrete_util::uuid_hash >;

  template <typename T>
  inline const T* const first_set_with_name(const std::vector<T>& obj_of_interest,
					    const std::string& tool_name) {
    for(typename std::vector< T >::const_iterator it = obj_of_interest.begin();
	it != obj_of_interest.end(); ++it) {
      if(it->metadata.tool.find(tool_name) != std::string::npos) return &(*it);  
    }
    return NULL;
  };
  template <typename T>
  inline const T* const first_set_with_name(const std::vector<T>& obj_of_interest,
					    const std::string& tool_name,
					    const std::string& type) {
    for(typename std::vector< T >::const_iterator it = obj_of_interest.begin();
	it != obj_of_interest.end(); ++it) {
      if(it->metadata.tool.find(tool_name) != std::string::npos &&
	 it->taggingType == type) return &(*it);
    }
    return NULL;
  };
  template <typename T>
  inline const T* const first_set(const std::vector<T>& obj_of_interest,
				  const std::string& type) {
    for(typename std::vector< T >::const_iterator it = obj_of_interest.begin();
	it != obj_of_interest.end(); ++it) {
      if(it->taggingType == type) return &(*it);
    }
    return NULL;
  };

  const concrete::TokenTagging* const first_pos_tagging(const concrete::Tokenization& tokenization,
							const std::string& tool_name);
  const concrete::TokenTagging* const first_ner_tagging(const concrete::Tokenization& tokenization,
							const std::string& tool_name);
  const concrete::TokenTagging* const first_lemma_tagging(const concrete::Tokenization& tokenization,
							  const std::string& tool_name);
  const concrete::DependencyParse* const first_dependency_parse(const concrete::Tokenization& tokenization,
								const std::string& tool_name);
  const concrete::EntitySet* const first_entity_set(const concrete::Communication& comm,
						    const std::string& tool_name);

  const uuid_map<concrete::EntityMention> mention_id_to_mention(const concrete::Communication& comm, const std::string& tool_name);
  const uuid_map<concrete::Tokenization> mention_id_to_tokenization(const concrete::Communication& comm, const std::string& tool_name);
  const uuid_map< std::list< concrete::SituationMention > > tokenization_id_to_situation_mention(const concrete::Communication& comm, const std::string& tool_name);
  uuid_map<concrete::Tokenization> tokenization_id_to_tokenization(const concrete::Communication& comm);

}

namespace ferrum {
  int get_num_mentions(const concrete::Communication& comm);
}

#include "ferrum/concrete.tcc"

#endif
