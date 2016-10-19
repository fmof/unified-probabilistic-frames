#ifndef FERRUM_LIBNAR_CRTLDA_CONCRETE_H_
#define FERRUM_LIBNAR_CRTLDA_CONCRETE_H_

#include "ferrum/concrete.hpp"
#include "ferrum/crtlda.hpp"
#include "ferrum/data_util.hpp"
#include "ferrum/tar.hpp"

#include "concrete_util/io.h"

#include <fcntl.h>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>
#include <unordered_map>

namespace ferrum {
  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol
    >
  inline void set_archived_corpus
  (
   const std::string& fpath,
   int N,
   unsigned int num_dep_hops,
   const ferrum::Toolnames& tools,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab,
   CorpusType* corpus
   ) {
    int num_comms = 0;
    ferrum::CompressedTar archive(fpath, ferrum::CompressedTarOperation::READ);
    for(ferrum::CompressedTar::value_type tar_entry : archive) {
      std::string entry_name(archive.name(tar_entry));
      TRACE << "Processing " << entry_name;
      concrete::Communication comm;
      concrete::util::ConcreteResultStruct crs =
	concrete::util::load_from_archive(archive, comm, ferrum::thrift::FromBuffer<Protocol>(), tar_entry);
      switch(crs.crr) {
      case concrete::util::ConcreteReadResult::OK:
	{
	  //comm.__set_id(entry_name);
	  if( corpus->has_document(comm.id) ) continue;
	  typename CorpusType::DocType doc = Pruner::make(comm, num_dep_hops, &gov_vocab, &rel_vocab, tools);
	  corpus->add_document(doc);
	  ++num_comms;
	}
	break;
      case concrete::util::ConcreteReadResult::E_BAD_SIZE:
	WARN << "Unable to read anything for communication " << entry_name;
	break;
      case concrete::util::ConcreteReadResult::E_OOM:
	ERROR << "Unable to malloc size " << crs.size << " while reading " << entry_name;
	throw std::bad_alloc();
      default:
	ERROR << "Unknown error";
	break;
      }
      if(num_comms == N) {
	INFO << "Loaded " << N << " communications, as requested";
	break;
      }
    } // end for loop
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol
    >
  inline void set_archived_corpus
  (
   const std::string& fpath,
   unsigned int num_dep_hops,
   const ferrum::Toolnames& tools,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab,
   CorpusType* corpus
   ) {
    int num_comms = 0;
    ferrum::CompressedTar archive(fpath, ferrum::CompressedTarOperation::READ);
    for(ferrum::CompressedTar::value_type tar_entry : archive) {
      std::string entry_name(archive.name(tar_entry));
      TRACE << "Processing " << entry_name;
      concrete::Communication comm;
      concrete::util::ConcreteResultStruct crs =
	concrete::util::load_from_archive(archive, comm, ferrum::thrift::FromBuffer<Protocol>(), tar_entry);
      switch(crs.crr) {
      case concrete::util::ConcreteReadResult::OK:
	{
	  //comm.__set_id(entry_name);
	  if( corpus->has_document(comm.id) ) continue;
	  typename CorpusType::DocType doc = Pruner::make(comm, num_dep_hops, &gov_vocab, &rel_vocab, tools);
	  corpus->add_document(doc);
	  ++num_comms;
	}
	break;
      case concrete::util::ConcreteReadResult::E_BAD_SIZE:
	WARN << "Unable to read anything for communication " << entry_name;
	break;
      case concrete::util::ConcreteReadResult::E_OOM:
	ERROR << "Unable to malloc size " << crs.size << " while reading " << entry_name;
	throw std::bad_alloc();
      default:
	ERROR << "Unknown error";
	break;
      }
    }
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename... MakerArgs
    >
  inline void set_archived_corpus
  (
   const std::string& fpath,
   int N,
   CorpusType* corpus,
   MakerArgs ... maker_args
   ) {
    int num_comms = 0;
    ferrum::CompressedTar archive(fpath, ferrum::CompressedTarOperation::READ);
    for(ferrum::CompressedTar::value_type tar_entry : archive) {
      std::string entry_name(archive.name(tar_entry));
      TRACE << "Processing " << entry_name;
      concrete::Communication comm;
      concrete::util::ConcreteResultStruct crs =
  	concrete::util::load_from_archive(archive, comm, ferrum::thrift::FromBuffer<Protocol>(), tar_entry);
      switch(crs.crr) {
      case concrete::util::ConcreteReadResult::OK:
  	{
  	  //comm.__set_id(entry_name);
  	  if( corpus->has_document(comm.id) ) continue;
  	  typename CorpusType::DocType doc = Pruner::make(comm, maker_args...);
  	  corpus->add_document(doc);
  	  ++num_comms;
  	}
  	break;
      case concrete::util::ConcreteReadResult::E_BAD_SIZE:
  	WARN << "Unable to read anything for communication " << entry_name;
  	break;
      case concrete::util::ConcreteReadResult::E_OOM:
  	ERROR << "Unable to malloc size " << crs.size << " while reading " << entry_name;
  	throw std::bad_alloc();
      default:
  	ERROR << "Unknown error";
  	break;
      }
      if(N == num_comms) {
  	INFO << "Loaded " << N << " communications, as requested";
  	break;
      }
    }
  }
  
  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename... MakerArgs
    >
  inline void set_archived_corpus
  (
   const std::string& fpath,
   CorpusType* corpus,
   MakerArgs ... maker_args
   ) {
    int num_comms = 0;
    ferrum::CompressedTar archive(fpath, ferrum::CompressedTarOperation::READ);
    for(ferrum::CompressedTar::value_type tar_entry : archive) {
      std::string entry_name(archive.name(tar_entry));
      TRACE << "Processing " << entry_name;
      concrete::Communication comm;
      concrete::util::ConcreteResultStruct crs =
	concrete::util::load_from_archive(archive, comm, ferrum::thrift::FromBuffer<Protocol>(), tar_entry);
      switch(crs.crr) {
      case concrete::util::ConcreteReadResult::OK:
	{
	  //comm.__set_id(entry_name);
	  if( corpus->has_document(comm.id) ) continue;
	  typename CorpusType::DocType doc = Pruner::make(comm, maker_args...);
	  corpus->add_document(doc);
	  ++num_comms;
	}
	break;
      case concrete::util::ConcreteReadResult::E_BAD_SIZE:
	WARN << "Unable to read anything for communication " << entry_name;
	break;
      case concrete::util::ConcreteReadResult::E_OOM:
	ERROR << "Unable to malloc size " << crs.size << " while reading " << entry_name;
	throw std::bad_alloc();
      default:
	ERROR << "Unknown error";
	break;
      }
    }
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol
    >
  inline CorpusType* get_archived_corpus
  (
   const std::string& name,
   const std::string& fpath,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab
   ) {
    CorpusType* corpus = new CorpusType(name);
    ferrum::Toolnames tools;
    set_archived_corpus<CorpusType, Pruner, Protocol>(fpath, 1, tools, gov_vocab, rel_vocab, corpus);
    return corpus;
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename ...Args
    >
  inline CorpusType* get_archived_corpus
  (
   const std::string& fpath,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab,
   Args... corpus_args
   ) {
    CorpusType* corpus = new CorpusType(corpus_args...);
    ferrum::Toolnames tools;
    set_archived_corpus<CorpusType, Pruner, Protocol>(fpath, 1, tools, gov_vocab, rel_vocab, corpus);
    return corpus;
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename ...Args
    >
  inline CorpusType* get_archived_corpus
  (
   const std::string& fpath,
   unsigned int num_dep_hops,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab,
   Args... corpus_args
   ) {
    CorpusType* corpus = new CorpusType(corpus_args...);
    ferrum::Toolnames tools;
    set_archived_corpus<CorpusType, Pruner, Protocol>(fpath, num_dep_hops, tools, gov_vocab, rel_vocab, corpus);
    return corpus;
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename ...Args
    >
  inline CorpusType* get_archived_corpus
  (
   const std::string& fpath,
   unsigned int num_dep_hops,
   ferrum::Vocabulary<std::string>& gov_vocab,
   ferrum::Vocabulary<std::string>& rel_vocab,
   const ferrum::Toolnames& tools,
   Args... corpus_args
   ) {
    CorpusType* corpus = new CorpusType(corpus_args...);
    set_archived_corpus<CorpusType, Pruner, Protocol>(fpath, num_dep_hops, tools, gov_vocab, rel_vocab, corpus);
    return corpus;
  }

  template <
    typename CorpusType,
    typename Pruner,
    typename Protocol = concrete::util::TCompactProtocol,
    typename ...Args
    >
  inline CorpusType* get_archived_corpus
  (
   const std::string& fpath,
   //std::tuple< MArgs... > maker_args,
   ferrum::Vocabulary<std::string>& vocab,
   const ferrum::Toolnames& tools,
   minsky::WordAnnotation::type wform,
   Args... corpus_args
   ) {
    CorpusType* corpus = new CorpusType(corpus_args...);
    set_archived_corpus<CorpusType, Pruner, Protocol>(fpath, corpus, &vocab, tools, wform);
    return corpus;
  }

}

#endif
