#include "ferrum/concrete.hpp"
#include "ferrum/thrift_protocol_defs.hpp"
#include <array>
#include <list>
#include <tuple>
#include <vector>
#include <unordered_map>

#include "ferrum/logging.hpp"

#include <boost/make_shared.hpp>
#include <chrono>
#include <cstring>
#include <errno.h>
#include <sys/stat.h>

namespace concrete {
  namespace util {
    /**
     * return true iff first is fully-contained in second
     */
    bool contained_in(const concrete::TextSpan& first,
		      const concrete::TextSpan& second) {
      return \
	(first.start >= second.start) && 
	(first.ending <= second.ending) &&
	(first.start <= second.ending) &&
	(first.ending >= second.start);
    }

    concrete::AnnotationMetadata make_metadata(const std::string& tool) {
      concrete::AnnotationMetadata am;
      am.__set_tool(tool);
      using namespace std::chrono;
      seconds secs = duration_cast< seconds >
	(
	 system_clock::now().time_since_epoch()
	 );
      long long lls = secs.count();
      am.__set_timestamp(lls);
      return am;
    }

    void get_binary_communication_sequence(const std::string& f_path_name, CommunicationSequence*& csp) {
      boost::filesystem::path f_path(f_path_name);
      if(boost::filesystem::is_directory(f_path)) {
	csp = new TBinaryDirectoryCommunicationSequence(f_path_name);
      } else {
	csp = new TBinaryGZipCommunicationSequence(f_path_name);
      }
    }

    template <typename P>
    ConcreteSmartWriter<P>::ConcreteSmartWriter(const std::string& base) : 
      ferrum::thrift::ThriftSmartWriter<P>(base),
      fd_(-1) {
    }

    template <typename P>
    P* ConcreteSmartWriter<P>::get(const std::string& suffix) {
      std::string f = this->next_file_name(suffix);
      if(fd_ >= 0) {
	INFO << "Closing previously opened file " << this->curr_file << " and opening a new one " << f;
	close_csw();
      }
      fd_ = open
	(
	 f.c_str(),
	 O_CREAT | O_WRONLY,
	 // set mode to rwxr--r--
	 (S_IRUSR | S_IWUSR | S_IXUSR) |
	 (S_IRGRP) |
	 (S_IROTH)
	 );
      if(fd_ < 0) {
	ERROR << "An error (" << errno << ", " <<  std::strerror(errno) << ") was encountered in trying to open " << f;
	throw 10;
      }
      this->curr_file = f;
      innerTransport_ = boost::make_shared<TFDTransport>(fd_);
      transport_ = boost::make_shared<TBufferedTransport>(innerTransport_);
      this->protocol_ = boost::make_shared<P>(transport_);
      transport_->open();
      return this->protocol_.get();
    }
    template <typename P>
    P* ConcreteSmartWriter<P>::get() {
      return this->get("");
    }
    template <typename P>
    int ConcreteSmartWriter<P>::close_csw() {
      int res = 0;
      if(fd_ >= 0) {
	transport_->close();
	res = close(fd_);
      }
      return res;
    }
    template <typename P>
    ConcreteSmartWriter<P>::~ConcreteSmartWriter() {
      close_csw();
    }

    template <typename P>
    void ConcreteSmartWriter<P>::_save() {
    }

    template <typename P>
    ConcreteSmartWriter<P>* ConcreteSmartWriter<P>::clone() const {
      ConcreteSmartWriter<P>* rtsw = new ConcreteSmartWriter<P>(*this);
      return rtsw;
    }

    ////////////////////////////////////////////////////////////////////

    // specialize
    template <>
    void StructIterator<concrete::Tokenization,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator >::ConcreteStructIteratorReturn::init_iters_() {
      auto& first = std::get<0>(iters_);
      first.c = par_->c_->sectionList.begin();
      first.e = par_->c_->sectionList.end();
      need_to_reset_[0] = false;
      //reset_inner_(index<1>());
      advance_(true, true);
    }
    template <>
    void StructIterator<concrete::Tokenization,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator >::ConcreteStructIteratorReturn::set_obj_ptr_() {
      auto& curr = std::get<1>(iters_);
      if(curr.c != curr.e) {
	obj_ptr_ = &(curr.c->tokenization);
      } else {
	obj_ptr_ = NULL;
      }
    }
    template <>
    template <>
    void StructIterator<concrete::Tokenization,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator >::ConcreteStructIteratorReturn::reset_inner_(index<1>) {
      constexpr size_t ct_which = 1;
      if( ! need_to_reset_[ct_which] ) {
      	return;
      } else if(ct_which > 0) {
      	auto& curr = std::get<ct_which>(iters_);
      	constexpr size_t prev_i = ct_which - 1;
      	auto& prev = std::get< prev_i >(iters_);
	if(prev.c != prev.e) {
	  curr.c = prev.c->sentenceList.begin();
	  curr.e = prev.c->sentenceList.end();
	  need_to_reset_[ct_which] = false;
	  while(curr.c != curr.e && !curr.c->__isset.tokenization) {
	    ++(curr.c);
	  }
	  if(curr.c != curr.e) {
	    set_obj_ptr_();
	  } else {
	    obj_ptr_ = NULL;
	    advance_(false, true);
	  }
	}
      }
    }

    // specialize: StructIterator<Communication>
    template <>
    void StructIterator<concrete::Communication, int>::ConcreteStructIteratorReturn::set_obj_ptr_() {
    }
    template <>
    void StructIterator<concrete::Communication, int>::ConcreteStructIteratorReturn::init_iters_() {
      auto& first = std::get<0>(iters_);
      first.c = 0;
      first.e = 1;
      need_to_reset_[0] = false;
      obj_ptr_ = par_->c_;
    }
    template <>
    template <size_t ct_which>
    void StructIterator<concrete::Communication, int>::ConcreteStructIteratorReturn::reset_inner_(index<ct_which>) {
    }
    template <>
    template <>
    void StructIterator<concrete::Communication, int >::ConcreteStructIteratorReturn::reset_inner_(index<0>) {
      set_obj_ptr_();
    }

    // specialize: StructIterator<Token>
    template <> template <>
    void StructIterator<concrete::Token,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator,
			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<0>) {
    }
    template <> template <>
    void StructIterator<concrete::Token,
    			std::vector<concrete::Section>::const_iterator,
    			std::vector<concrete::Sentence>::const_iterator,
    			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<1>) {
      constexpr size_t ct_which = 1;
      reset_inner_(index<0>());
      if( ! need_to_reset_[ct_which] ) {
      	return;
      } else {
      	auto& curr = std::get<ct_which>(iters_);
    	constexpr size_t prev_i = ct_which - 1;
    	auto& prev = std::get< prev_i >(iters_);
	while(true) {
	  // are there any sentences in this section?
	  // if not, advance the section counter, until we
	  // find sentences or we're out of sections
	  while(prev.c != prev.e &&
		( (curr.c = prev.c->sentenceList.begin()) ==
		  (curr.e = prev.c->sentenceList.end()) ) ) {
	    ++(prev.c);
	  }
	  // it could be we've exhausted the sections
	  if(prev.c == prev.e) {
	    obj_ptr_ = NULL;
	    return;
	  }
	  // otherwise, iterate through the sentences until we've found one
	  // with a set tokenization
	  while(curr.c != curr.e && !curr.c->__isset.tokenization) {
	    ++(curr.c);
	  }
	  if(curr.c != curr.e) {
	    need_to_reset_[ct_which] = false;
	    break;
	  } else {
	    // continue w/ the loop
	  }
	}
      }
    }
    template <>
    void StructIterator<concrete::Token,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator,
			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::set_obj_ptr_() {
      auto& curr = std::get<2>(iters_);
      if(curr.c != curr.e) {
	obj_ptr_ = &(*curr.c);
      } else {
	obj_ptr_ = NULL;
      }
    }
    template <> template <>
    void StructIterator<concrete::Token,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator,
			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<2>) {
      constexpr size_t ct_which = 2;
      reset_inner_(index<1>());
      if( ! need_to_reset_[ct_which] ) {
      	return;
      } else if( need_to_reset_[1] ) {
	return ;
      } else {
      	auto& curr = std::get<ct_which>(iters_);
	constexpr size_t prev_i = ct_which - 1;
	auto& prev = std::get< prev_i >(iters_);
	while( (prev.c != prev.e) &&
	       ((curr.c = prev.c->tokenization.tokenList.tokenList.begin()) ==
		(curr.e = prev.c->tokenization.tokenList.tokenList.end() ) ) ) {
	  ++(prev.c);
	}
	if(prev.c == prev.e) { // need to update the grandchild iterator
	  int gp_res = check<0>(false);
	  if(gp_res < 0) {
	    obj_ptr_ = NULL;
	    return ;
	  }
	  need_to_reset_[1] = true;
	  reset_inner_(index<2>());
	} else {
	  if(curr.c != curr.e) {
	    obj_ptr_ = &(*curr.c);
	    need_to_reset_[ct_which] = false;
	  } else {
	    obj_ptr_ = NULL;
	    advance_(false, true);
	  }
	}
      }
    }

    template <>
    void StructIterator<concrete::Token,
			std::vector<concrete::Section>::const_iterator,
			std::vector<concrete::Sentence>::const_iterator,
			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::init_iters_() {
      auto& first = std::get<0>(iters_);
      first.c = par_->c_->sectionList.begin();
      first.e = par_->c_->sectionList.end();
      need_to_reset_[0] = false;
      //reset_inner_(index<2>());
      advance_(true, true);
    }

    // // specialize: StructIterator<TaggedToken>
    // template <> template <>
    // void StructIterator<concrete::TaggedToken,
    // 			std::vector<concrete::Section>::const_iterator,
    // 			std::vector<concrete::Sentence>::const_iterator,
    // 			std::vector<concrete::TaggedToken>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<0>) {
    // }
    // template <> template <>
    // void StructIterator<concrete::TaggedToken,
    // 			std::vector<concrete::Section>::const_iterator,
    // 			std::vector<concrete::Sentence>::const_iterator,
    // 			std::vector<concrete::TaggedToken>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<1>) {
    //   constexpr size_t ct_which = 1;
    //   reset_inner_(index<0>());
    //   if( ! need_to_reset_[ct_which] ) {
    //   	return;
    //   } else {
    //   	auto& curr = std::get<ct_which>(iters_);
    // 	constexpr size_t prev_i = ct_which - 1;
    // 	auto& prev = std::get< prev_i >(iters_);
    // 	while(true) {
    // 	  // are there any sentences in this section?
    // 	  // if not, advance the section counter, until we
    // 	  // find sentences or we're out of sections
    // 	  while(prev.c != prev.e &&
    // 		( (curr.c = prev.c->sentenceList.begin()) ==
    // 		  (curr.e = prev.c->sentenceList.end()) ) ) {
    // 	    ++(prev.c);
    // 	  }
    // 	  // it could be we've exhausted the sections
    // 	  if(prev.c == prev.e) {
    // 	    obj_ptr_ = NULL;
    // 	    return;
    // 	  }
    // 	  // otherwise, iterate through the sentences until we've found one
    // 	  // with a set tokenization
    // 	  while(curr.c != curr.e && !curr.c->__isset.tokenization) {
    // 	    ++(curr.c);
    // 	  }
    // 	  if(curr.c != curr.e) {
    // 	    need_to_reset_[ct_which] = false;
    // 	    break;
    // 	  } else {
    // 	    // continue w/ the loop
    // 	  }
    // 	}
    //   }
    // }
    // template <>
    // void StructIterator<concrete::TaggedToken,
    // 			std::vector<concrete::Section>::const_iterator,
    // 			std::vector<concrete::Sentence>::const_iterator,
    // 			std::vector<concrete::TaggedToken>::const_iterator>::ConcreteStructIteratorReturn::set_obj_ptr_() {
    //   auto& curr = std::get<2>(iters_);
    //   if(curr.c != curr.e) {
    // 	obj_ptr_ = &(*curr.c);
    //   } else {
    // 	obj_ptr_ = NULL;
    //   }
    // }
    // template <> template <>
    // void StructIterator<concrete::TaggedToken,
    // 			std::vector<concrete::Section>::const_iterator,
    // 			std::vector<concrete::Sentence>::const_iterator,
    // 			std::vector<concrete::TaggedToken>::const_iterator>::ConcreteStructIteratorReturn::reset_inner_(index<2>) {
    //   constexpr size_t ct_which = 2;
    //   reset_inner_(index<1>());
    //   if( ! need_to_reset_[ct_which] ) {
    //   	return;
    //   } else if( need_to_reset_[1] ) {
    // 	return ;
    //   } else {
    //   	auto& curr = std::get<ct_which>(iters_);
    // 	constexpr size_t prev_i = ct_which - 1;
    // 	auto& prev = std::get< prev_i >(iters_);
    // 	while( (prev.c != prev.e) &&
    // 	       ((curr.c = prev.c->tokenization.tokenList.tokenList.begin()) ==
    // 		(curr.e = prev.c->tokenization.tokenList.tokenList.end() ) ) ) {
    // 	  ++(prev.c);
    // 	}
    // 	if(prev.c == prev.e) { // need to update the grandchild iterator
    // 	  int gp_res = check<0>(false);
    // 	  if(gp_res < 0) {
    // 	    obj_ptr_ = NULL;
    // 	    return ;
    // 	  }
    // 	  need_to_reset_[1] = true;
    // 	  reset_inner_(index<2>());
    // 	} else {
    // 	  if(curr.c != curr.e) {
    // 	    obj_ptr_ = &(*curr.c);
    // 	    need_to_reset_[ct_which] = false;
    // 	  } else {
    // 	    obj_ptr_ = NULL;
    // 	    advance_(false, true);
    // 	  }
    // 	}
    //   }
    // }

    // template <>
    // void StructIterator<concrete::Token,
    // 			std::vector<concrete::Section>::const_iterator,
    // 			std::vector<concrete::Sentence>::const_iterator,
    // 			std::vector<concrete::Token>::const_iterator>::ConcreteStructIteratorReturn::init_iters_() {
    //   auto& first = std::get<0>(iters_);
    //   first.c = par_->c_->sectionList.begin();
    //   first.e = par_->c_->sectionList.end();
    //   need_to_reset_[0] = false;
    //   //reset_inner_(index<2>());
    //   advance_(true, true);
    // }

  } // ends namespace util
} //ends namespace concrete

template class concrete::util::ConcreteSmartWriter<ferrum::thrift::TBinaryProtocol>;
template class concrete::util::ConcreteSmartWriter<ferrum::thrift::TCompactProtocol>;
template class concrete::util::ConcreteSmartWriter<ferrum::thrift::TJSONProtocol>;


namespace concrete_util {
  const size_t uuid_hash::operator()(const concrete::UUID& uuid) const {
    return std::hash<std::string>()(uuid.uuidString);
  }

  const concrete::TokenTagging* const first_pos_tagging(const concrete::Tokenization& tokenization,
							const std::string& tool_name) {
    return concrete_util::first_set_with_name<concrete::TokenTagging>(tokenization.tokenTaggingList, tool_name, "POS");
  }
  const concrete::TokenTagging* const first_ner_tagging(const concrete::Tokenization& tokenization,
							const std::string& tool_name) {
    return concrete_util::first_set_with_name<concrete::TokenTagging>(tokenization.tokenTaggingList, tool_name, "NER");
  }
  const concrete::TokenTagging* const first_lemma_tagging(const concrete::Tokenization& tokenization,
							  const std::string& tool_name) {
    return concrete_util::first_set_with_name<concrete::TokenTagging>(tokenization.tokenTaggingList, tool_name, "LEMMA");
  }

  const concrete::DependencyParse* const first_dependency_parse(const concrete::Tokenization& tokenization,
								const std::string& tool_name) {
    return concrete_util::first_set_with_name<concrete::DependencyParse>(tokenization.dependencyParseList, tool_name);
  }

  const concrete::EntitySet* const first_entity_set(const concrete::Communication& comm,
						    const std::string& tool_name) {
    return concrete_util::first_set_with_name<concrete::EntitySet>(comm.entitySetList, tool_name);
  }

  uuid_map<concrete::Tokenization> tokenization_id_to_tokenization(const concrete::Communication& comm) {
    uuid_map<concrete::Tokenization> tutt;
    int num_tok = 0;
    if(!comm.__isset.sectionList) return tutt;
    for(const concrete::Section& section : comm.sectionList) {
      if(!section.__isset.sentenceList) continue;
      //for(const concrete::Sentence& sentence : section.sentenceList) {
      for(std::vector<concrete::Sentence>::const_iterator sit = section.sentenceList.begin();
	  sit != section.sentenceList.end(); ++sit) {
	if(!sit->__isset.tokenization) continue;
	const concrete::Tokenization* tptr = &(sit->tokenization);
	//DEBUG << "Storing tokenization idx = " << num_tok << " with ID " << tptr->uuid.uuidString;
	tutt[tptr->uuid] = *tptr;
	//DEBUG << "checking: " << tutt[tptr->uuid].uuid.uuidString;
	//DEBUG << "...loaded";
	++num_tok;
      }
    }
    return tutt;
  }

  const uuid_map<concrete::EntityMention> mention_id_to_mention(const concrete::Communication& comm, const std::string& tool_name) {
    const concrete::EntityMentionSet* ems =
      concrete_util::first_set_with_name<concrete::EntityMentionSet>(comm.entityMentionSetList,
								     tool_name);
    // select the theory corresponding to toolname
    uuid_map<concrete::EntityMention> mitm;
    if(!ems) return mitm;
    for(concrete::EntityMention em : ems->mentionList) {
      mitm[em.uuid] = em;
    }
    return mitm;
  }

  const uuid_map< std::list< concrete::SituationMention > > tokenization_id_to_situation_mention(const concrete::Communication& comm, const std::string& tool_name) {
    const concrete::SituationMentionSet* sms =
      concrete_util::first_set_with_name<concrete::SituationMentionSet>(comm.situationMentionSetList,
									tool_name);
    uuid_map< std::list< concrete::SituationMention> > tism;
    if(!sms) {
      WARN << "Did not find any SituationMentionSets in communication " << comm.id << " with name containing " << tool_name;
      return tism;
    }
    for(const concrete::SituationMention sm : sms->mentionList) {
      // if(tism.find(sm.tokens.tokenizationId) == tism.end()) {
      // 	tism[sm.tokens.tokenizationId] = std::list<concrete::SituationMention>();
      // }
      tism[sm.tokens.tokenizationId].push_back(sm);
    }
    return tism;
  }

  const uuid_map<concrete::Tokenization> mention_id_to_tokenization(const concrete::Communication& comm, const std::string& tool_name) {
    const concrete::EntityMentionSet* ems =
      concrete_util::first_set_with_name<concrete::EntityMentionSet>(comm.entityMentionSetList,
								     tool_name);
    uuid_map<concrete::Tokenization> mitt;
    if(!ems) {
      WARN << "Did not find any EntityMentionSets in communication " << comm.id << " with name containing " << tool_name;
      return mitt;
    }
    uuid_map<concrete::Tokenization> tid_to_tok =
      concrete_util::tokenization_id_to_tokenization(comm);
    for(const concrete::EntityMention em : ems->mentionList) {
      concrete::Tokenization tp = tid_to_tok.at(em.tokens.tokenizationId);
      mitt[em.uuid] = tp;
    }
    return mitt;
  }
}

namespace ferrum {
  int get_num_mentions(const concrete::Communication& comm) {
    int num = 0;
    for(const concrete::EntityMentionSet& ems : comm.entityMentionSetList) {
      num += ems.mentionList.size();
    }
    return num;
  }
}
