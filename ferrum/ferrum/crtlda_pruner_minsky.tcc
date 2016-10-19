#ifndef FERRUM_LIBNAR_CRTLDA_MINSKY_PRUNER_TCC_
#define FERRUM_LIBNAR_CRTLDA_MINSKY_PRUNER_TCC_

#include "ferrum/util.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/crtlda_pruner_base.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/tm_pruner.hpp"

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include <time.h>

// for pair
#include <map>
#include <utility>
#include <unordered_set>
#include <string>
#include <vector>

namespace ferrum {
  template <typename Subpruner>
  void MinskyPruner<Subpruner>::make_edoc(minsky::EDoc* doc) {
    doc->__set_id(communication.id);
    const concrete::EntitySet* const esp =
      concrete_util::first_entity_set(communication,
				      tools_.entity_mention_tool);
    if(! esp) return;
    const std::vector< concrete::Entity >& c_entity_set = esp->entityList;
    concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention =
      concrete_util::mention_id_to_mention
      (
       communication,
       tools_.entity_mention_tool
       );
    size_t total_num_mentions = 0;
    size_t num_entities = 0;
    typedef std::pair<std::string, concrete::UUID> str_pair;
    std::map< std::pair<std::string, int>, str_pair > filter_map;
    for(const concrete::Entity& conc_entity : c_entity_set) {
      for(const concrete::UUID& mention_uuid : conc_entity.mentionIdList) {
	const concrete::EntityMention& conc_entity_mention =
	  mention_id_to_mention[mention_uuid];
	const concrete::Tokenization& mention_tokenization =
	  this->mention_tokenization(mention_uuid);
	const int anchor_index = conc_entity_mention.tokens.anchorTokenIndex;
	std::pair<std::string, int> key(mention_tokenization.uuid.uuidString, anchor_index);
	// filter_map[key] = str_pair(conc_entity.uuid.uuidString,
	// 			   mention_uuid);
	auto it = filter_map.find(key);
	if(it == filter_map.end()) {
	  filter_map[key] = str_pair(conc_entity.uuid.uuidString,
				     mention_uuid);
	} else {
	  const str_pair& previous = it->second;
	  const concrete::EntityMention& prev_mention = mention_id_to_mention[previous.second];
	  if(conc_entity_mention.tokens.__isset.textSpan && prev_mention.tokens.__isset.textSpan) {
	    const concrete::TextSpan& curr_text_span = conc_entity_mention.tokens.textSpan;
	    const concrete::TextSpan& prev_text_span = prev_mention.tokens.textSpan;
	    // now if we've found a larger span, then take that
	    if(concrete::util::contained_in(prev_text_span, curr_text_span) ) {
	      INFO << "Replacing concrete::EntityMention " << prev_mention.uuid.uuidString << " with " << mention_uuid.uuidString << \
		": (" << prev_text_span.start << ", " << prev_text_span.ending << ") is contained in " << \
		"(" << curr_text_span.start << ", " << curr_text_span.ending << ")";
	      filter_map[key] =
		str_pair(conc_entity.uuid.uuidString,
			 mention_uuid);
	    }
	  } else {
	    const size_t curr_span = conc_entity_mention.tokens.tokenIndexList.size();
	    const size_t prev_span = prev_mention.tokens.tokenIndexList.size();
	    bool replace_span_overlap = ferrum::contained_in
	      (ferrum::span(prev_mention.tokens.tokenIndexList),
	       ferrum::span(conc_entity_mention.tokens.tokenIndexList));
	    if(replace_span_overlap ||
	       (curr_span > prev_span)) {
	      filter_map[key] =
		str_pair(conc_entity.uuid.uuidString,
			 mention_uuid);
	      INFO << "Replacing concrete::EntityMention " << prev_mention.uuid.uuidString << " with " << mention_uuid.uuidString << \
		   " due to " << (replace_span_overlap ? " overlapping spans" : "longer, non-subsumed spans");
	    }
	  }
	}
      }
    } // end filter_map construction loop
    for(const concrete::Entity& conc_entity : c_entity_set) {
      minsky::Entity entity;
      //entity.canonical_name(conc_entity.canonicalName);
      entity.__set_id(conc_entity.uuid.uuidString);
      size_t num_mentions = 0;
      for(const concrete::UUID& mention_uuid : conc_entity.mentionIdList) {
	const concrete::EntityMention& conc_entity_mention =
	  mention_id_to_mention[mention_uuid];
	const concrete::Tokenization& mention_tokenization =
	  this->mention_tokenization(mention_uuid);
	const int anchor_index = conc_entity_mention.tokens.anchorTokenIndex;
	std::pair<std::string, int> key(mention_tokenization.uuid.uuidString, anchor_index);
	const auto it = filter_map.find(key);
	if(it == filter_map.end()) {
	  ERROR << "After constructing the filter map, we can't find " << key.first << ", " << key.second;
	  throw 19;
	}
	if(it->second.first != conc_entity.uuid.uuidString ||
	   it->second.second.uuidString != mention_uuid.uuidString) {
	  continue;
	}
	std::vector< minsky::Mention > created_mentions =
	  this->prune(conc_entity_mention);
	if(created_mentions.size() > 0) {
	  DEBUG << "Adding " << created_mentions.size() << " mentions to " << entity.mentions.size() << " existing ones";
	  auto back_end = entity.mentions.end();
	  entity.mentions.insert(back_end, created_mentions.begin(), created_mentions.end());
	  entity.__isset.mentions = true;
	  num_mentions += created_mentions.size();
	}
      }
      if(num_mentions > 0) {
	DEBUG << "Entity " << entity.id << " has " << entity.mentions.size() << " mentions (doc " << doc->id << ")";
	doc->entities.push_back(entity);
	doc->__isset.entities = true;
	total_num_mentions += num_mentions;
	num_entities++;
      }
    }
    INFO << "Document " << doc->id << " has " << num_entities << " entities and " << total_num_mentions << " mentions";
  }
  
  template <typename Subpruner>
  MinskyPruner<Subpruner>::MinskyPruner
  (
   const concrete::Communication& comm,
   const ferrum::Toolnames tools
   ) :
    EntityMentionPruner< minsky::Mention >(comm), tools_(tools) {
  }
  template <typename Subpruner>
  MinskyPruner<Subpruner>::~MinskyPruner() {
  }

  template <typename Subpruner>
  template <typename ...ArgTypes>
  minsky::EDoc MinskyPruner<Subpruner>::make
  (
   const concrete::Communication& comm,
   unsigned int num_dep_hops,
   ArgTypes... args) {
    Subpruner sp(comm, args...);
    sp.num_dep_hops(num_dep_hops);
    TRACE << comm.id;
    minsky::EDoc doc;
    sp.make_edoc(&doc);
    return doc;
  }

  template <typename Subpruner>
  template <typename ...ArgTypes>
  minsky::EDoc MinskyPruner<Subpruner>::make(const concrete::Communication& comm,
					     ArgTypes... args) {
    return make<ArgTypes...>(comm, 1, args...);
  }
  template <typename Subpruner>
  void MinskyPruner<Subpruner>::num_dep_hops(unsigned int ndh) {
  }

  ////////////////////////////////////////////////////////////////////////

  template <typename Subpruner, typename CT>
  void MinskySentencePruner<Subpruner, CT>::make_simple_doc(minsky::SimpleDoc* doc) {
    doc->__set_id(this->communication.id);
    int num_clauses = 0;
    for(const CT& obj : concrete::util::StructIterator<CT>(this->communication)) {
      minsky::WordsClause cc = this->clause_create(obj);
      if(! minsky::empty(cc)) {
    	doc->sentences.push_back(cc);
    	doc->__isset.sentences = true;
    	num_clauses++;
      }
    }
    INFO << "Document " << doc->id << " has " << num_clauses << " clauses";
  }
  
  template <typename Subpruner, typename CT>
  MinskySentencePruner<Subpruner, CT>::MinskySentencePruner
  (
   const concrete::Communication& comm,
   const ferrum::Toolnames tools,
   minsky::WordAnnotation::type wannot
   ) :
    BOWPruner< minsky::WordsClause, CT >(comm), tools_(tools), annot_(wannot) {
  }
  template <typename Subpruner, typename CT>
  MinskySentencePruner<Subpruner, CT>::~MinskySentencePruner() {
  }

  template <typename Subpruner, typename CT>
  template <typename ...ArgTypes>
  minsky::SimpleDoc MinskySentencePruner<Subpruner, CT>::make
  (
   const concrete::Communication& comm,
   ArgTypes... args
   ) {
    //return make<ArgTypes...>(comm, 1, args...);
    Subpruner sp(comm, args...);
    //sp.num_dep_hops(num_dep_hops);
    TRACE << comm.id;
    minsky::SimpleDoc doc;
    sp.make_simple_doc(&doc);
    return doc;
  }
}

#endif
