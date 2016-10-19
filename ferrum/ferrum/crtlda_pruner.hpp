#ifndef FERRUM_LIBNAR_CRTLDA_PRUNER_H_
#define FERRUM_LIBNAR_CRTLDA_PRUNER_H_

#include "ferrum/dmc.hpp"
#include "ferrum/util.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_pruner_base.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <fstream>
#include <iostream>
#include "stdlib.h"
#include <time.h>

// for pair
#include "map"
#include <utility>
#include <unordered_set>
#include <string>
#include <vector>

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>

namespace ferrum {
  template <typename G, typename R>
  class VerbGovernedPruner : public EntityMentionPruner<Mention<G, R> > {
  private:
    const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization;
    std::string dep_parse_toolname = "col-ccproc-deps";
    typedef Mention<G,R> M;
  public:
    VerbGovernedPruner< G, R >(const concrete::Communication& comm) : EntityMentionPruner< Mention< G, R> >(comm), mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")) {
    }
    virtual const R make_relation(const std::string& concrete_relation, G gov_view) const {
      std::string str_gov_view(gov_view);
      std::string ret_str;
      ret_str = concrete_relation + "-" + str_gov_view;
      return ret_str;
    }
    virtual const G make_gov_view(const std::string& lemma) const {
      return lemma;
    }
    virtual const G make_head_view(const std::string& lemma) const {
      return lemma;
    }
    virtual std::vector<M> prune(const concrete::EntityMention& conc_mention) const {
      const concrete::Tokenization mention_tokenization = mention_id_to_tokenization.at(conc_mention.uuid);
      const int anchor_index = conc_mention.tokens.anchorTokenIndex;
      std::vector<M> created_mentions;
      if(anchor_index < 0) {
	INFO << "Anchor index not set in Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
	return created_mentions;
      }
      //get the dependency parse
      const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(mention_tokenization, dep_parse_toolname);
      if(dep_parse == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " does not have any dependency parses with name containing " << dep_parse_toolname;
	return created_mentions;
      }
      // for every dependency, check to see if dependency.dep == anchor_index
      for(const concrete::Dependency dependency : dep_parse->dependencyList) {
	if(dependency.dep != anchor_index || dependency.gov < 0) continue;
	const concrete::TokenTagging * const posTagging = concrete_util::first_pos_tagging(mention_tokenization, "Stanford");
	if(posTagging == NULL) {
	  WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	  continue;
	}
	// TODO: this really should be "generalizing"
	if(posTagging->taggedTokenList[dependency.gov].tag[0] != 'V') continue;
	const concrete::TokenTagging * const lemmaTagging = concrete_util::first_lemma_tagging(mention_tokenization, "Stanford");
	if(lemmaTagging == NULL) {
	  WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	  continue;
	}
	M mention;
	mention.latent(false);
	const std::string gov_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	AnnotatedToken< G > gov(gov_lemma, mention_tokenization.tokenList.tokenList[dependency.gov].text, posTagging->taggedTokenList[dependency.gov].tag);
	gov.view( this->make_gov_view(gov_lemma));
	const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	AnnotatedToken< G > head(lemmaTagging->taggedTokenList[dependency.dep].tag, mention_tokenization.tokenList.tokenList[dependency.dep].text, posTagging->taggedTokenList[dependency.dep].tag);
	head.view( this->make_head_view(head_lemma) );
	mention.gov(gov);
	mention.head(head);
	R rel = this->make_relation(dependency.edgeType, gov.view());
	mention.rel( rel );
	created_mentions.push_back(mention);
      }
      return created_mentions;
    }
  };

  template <typename G, typename R>
  class SituationGovernedPruner : public EntityMentionPruner<Mention<G, R> > {
  protected:
    typedef Vocabulary<G> GVocab;
    typedef Vocabulary<R> RVocab;
  private:
    // Record entity mention IDs to tokenizations
    const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization;
    // Record tokenizations to Situation mentions
    const concrete_util::uuid_map< std::list< concrete::SituationMention > > tok_id_to_sm;
    concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention;
     
    std::string dep_parse_toolname = "col-ccproc-deps";
    GVocab* gvocab_;
    RVocab* rvocab_;
    Vocabulary<std::string>* g_lemma_vocab_;
    Vocabulary<std::string>* r_lemma_vocab_;
    typedef Mention<G,R> M;
  protected:
    bool use_lexical_;
    bool add_to_vocabs_;
    StopWordList gov_stopwords_;
    bool use_gov_stopwords_;
  public:
    SituationGovernedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc, 
				    Vocabulary<std::string>* g_lemma_vocab,
				    Vocabulary<std::string>* r_lemma_vocab,
				    bool use_lex) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(g_lemma_vocab), r_lemma_vocab_(r_lemma_vocab),
      use_lexical_(use_lex), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }
    SituationGovernedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(false), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }
    SituationGovernedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc, bool use_lex) : 
    EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(use_lex), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }

    SituationGovernedPruner< G, R >(const concrete::Communication& comm) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(false), add_to_vocabs_(false), use_gov_stopwords_(false) {
    }

    void use_gov_stopwords(bool b) {
      use_gov_stopwords_ = b;
    }
    void gov_stopwords(const StopWordList& list_) {
      gov_stopwords_ = list_;
    }

    inline void gov_lemma_vocab(Vocabulary<std::string>* gv) {
      g_lemma_vocab_ = gv;
    }
    inline void rel_lemma_vocab(Vocabulary<std::string>* rv) {
      r_lemma_vocab_ = rv;
    }

    inline void use_lexical(bool b) {
      use_lexical_ = b;
    }
    bool use_lexical() {
      return use_lexical_;
    }
    bool use_lexical() const  {
      return use_lexical_;
    }
    virtual const R make_relation(const std::string rel, const std::string gov) const {
      std::string rel_ret = rel;
      rel_ret += "-" + gov;
      return rel_ret;
    }
    virtual const R make_relation(const std::string& concrete_relation, const std::string& dependency_rel,
				  G gov_view, const std::string& gov_lemma) const {
      return use_lexical_ ? make_relation(std::string(dependency_rel), std::string(gov_lemma)) : 
	make_relation(std::string(concrete_relation), std::string(gov_view));
    }
    virtual const G make_gov_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      return std::string(use_lexical_ ? token.lemma() : view_str);
    }
    virtual const G make_head_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      return std::string(use_lexical_ ? token.lemma() : view_str);
    }
    virtual std::vector<M> prune(const concrete::EntityMention& conc_mention) const {
      const concrete::Tokenization mention_tokenization = mention_id_to_tokenization.at(conc_mention.uuid);
      const int anchor_index = conc_mention.tokens.anchorTokenIndex;
      std::vector<M> created_mentions;
      if(tok_id_to_sm.find(mention_tokenization.uuid) == tok_id_to_sm.end()) {
	DEBUG << "No situation mentions found for Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
	return created_mentions;
      }
      if(anchor_index < 0) {
	DEBUG << "Anchor index not set in Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
	// we should probably log a warning here
	return created_mentions;
      }
      //get the dependency parse
      const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(mention_tokenization, dep_parse_toolname);
      if(dep_parse == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " does not have any dependency parses with name containing " << dep_parse_toolname;
	return created_mentions;
      }
      //get the SituationMention corresponding to this mention's tokenization
      std::list<concrete::SituationMention> sit_ment_list = tok_id_to_sm.at(mention_tokenization.uuid);
      const concrete::TokenTagging * const posTagging = concrete_util::first_pos_tagging(mention_tokenization, "Stanford");
      if(posTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	return created_mentions;
      }
      const concrete::TokenTagging * const lemmaTagging = concrete_util::first_lemma_tagging(mention_tokenization, "Stanford");
      if(lemmaTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no LEMMA tagging";
	return created_mentions;
      }
      // for every dependency, check to see if dependency.dep == anchor_index
      for(const concrete::Dependency dependency : dep_parse->dependencyList) {
	if(dependency.dep != anchor_index || dependency.gov < 0) continue;
	// TODO: this really should be "generalizing"
	if(posTagging->taggedTokenList[dependency.gov].tag[0] != 'V') continue;
	if(use_gov_stopwords_ && 
	   gov_stopwords_.contains(lemmaTagging->taggedTokenList[dependency.gov].tag)) {
	  continue;
	}
	// now verify that the governor corresponds to a situation mention trigger
	for(const concrete::SituationMention sit_ment : sit_ment_list) {
	  if(situation_mention_violations(sit_ment) || 
	     !concrete::util::index_overlap_mention<concrete::SituationMention>(sit_ment, dependency.gov)) continue;
	  // now, make sure that the dependency token is covered by one of the arcs
	  int mention_arg_idx = -1;
	  int mention_counter = -1;
	  for(const concrete::MentionArgument mention_argument : sit_ment.argumentList) {
	    ++mention_counter;
	    if(!mention_argument.__isset.role) {
	      WARN << "Skipping mention argument " << mention_counter << " due to missing `role` field";
	      continue;
	    }
	    if(!mention_argument.__isset.entityMentionId) {
	      WARN << "Skipping mention argument " << mention_counter << " due to missing `entityMentionId` field (this object only works with EntityMention arguments)";
	      continue;
	    }
	    // const concrete::Tokenization arg_mention_tokenization = 
	    //   sm_arg_mention_id_to_tokenization.at(mention_argument.entityMentionId);
	    const concrete::EntityMention arg_ent_mention =
	      mention_id_to_mention.at(mention_argument.entityMentionId);
	    if(!concrete::util::index_overlap_mention<concrete::EntityMention>(arg_ent_mention, dependency.dep)) continue;
	    mention_arg_idx = mention_counter;
	    break;
	  }
	  if(mention_arg_idx < 0) continue;
	  const concrete::MentionArgument rel_mention_arg = sit_ment.argumentList[mention_arg_idx];
	  M mention;
	  const std::string g_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	  AnnotatedToken< G > gov(g_lemma,
				  mention_tokenization.tokenList.tokenList[dependency.gov].text,
				  posTagging->taggedTokenList[dependency.gov].tag);
	  // The gov view is 
	  const std::string frame_trigger = sit_ment.situationKind;
	  const G gov_view = this->make_gov_view(gov, frame_trigger);
	  gov.view( gov_view );
	  mention.latent(false);
	  if(add_to_vocabs_) {
	    gvocab_->make_word(gov_view);
	    if(g_lemma_vocab_ != NULL) {
	      g_lemma_vocab_->make_word(g_lemma);
	    }
	  }
	  const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	  AnnotatedToken< G > head(head_lemma,
				   mention_tokenization.tokenList.tokenList[dependency.dep].text,
				   posTagging->taggedTokenList[dependency.dep].tag);
	  head.view( this->make_head_view(head, head_lemma) );
	  mention.gov(gov);
	  mention.head(head);
	  const R rel = this->make_relation(rel_mention_arg.role, dependency.edgeType, gov.view(), g_lemma);
	  const std::string r_lemma = this->make_relation(dependency.edgeType, g_lemma);
	  mention.rel_str( r_lemma );
	  mention.rel( rel );
	  if(add_to_vocabs_) {
	    rvocab_->make_word(rel);
	    if(r_lemma_vocab_ != NULL) {
	      r_lemma_vocab_->make_word(r_lemma);
	    }
	  }
	  mention.id(conc_mention.uuid.uuidString);
	  created_mentions.push_back(mention);
	}
      }
      return created_mentions;
    }

    /**
     * Return the number of assumption violations.
     * Any return value > 0 means we will skip processing this mention.
     */
    int situation_mention_violations(const concrete::SituationMention& sm) const {
      int num = 0;
      if(! sm.__isset.tokens) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to missing 'tokens' field";
	++num;
      }
      if(! sm.__isset.situationKind) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to missing 'situationKind' field";
	++num;
      }
      if(!sm.argumentList.size()) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to empty 'mentionArgument' field";
	++num;
      }
      return num;
    }
  };

  template <typename G, typename R>
  class SituationGovernedPrunerMaskLemma : public SituationGovernedPruner<G, R> {
  private:
    // add masking map to here
    WordMapper word_mapper_;
    typedef Vocabulary<G> GVocab;
    typedef Vocabulary<R> RVocab;
  public:
    SituationGovernedPrunerMaskLemma< G, R >(const concrete::Communication& comm,
					     GVocab *gvoc, RVocab* rvoc, 
					     Vocabulary<std::string>* g_lemma_vocab,
					     Vocabulary<std::string>* r_lemma_vocab,
					     bool use_lex,
					     const WordMapper& word_mapper) : SituationGovernedPruner< G, R >(comm, gvoc, rvoc, g_lemma_vocab, r_lemma_vocab, use_lex), word_mapper_(word_mapper) {
    }
    SituationGovernedPrunerMaskLemma< G, R >(const concrete::Communication& comm,
					     GVocab *gvoc, RVocab* rvoc, bool use_lex,
					     const WordMapper& word_mapper) : SituationGovernedPruner< G, R >(comm, gvoc, rvoc, use_lex), word_mapper_(word_mapper) {
    }
    SituationGovernedPrunerMaskLemma< G, R >(const concrete::Communication& comm,
					     const WordMapper& word_mapper) : SituationGovernedPruner< G, R >(comm), word_mapper_(word_mapper) {
    }

    /**
     * For this, we will always use the gov view, since that contains the mask transformaton
     */
    virtual const R make_relation(const std::string& concrete_relation, const std::string& dependency_rel,
    				  G gov_view, const std::string& gov_lemma) const {
      //return this->use_lexical() ? 
      return SituationGovernedPruner<G,R>::make_relation(std::string(dependency_rel), std::string(gov_view));
      //SituationGovernedPruner<G,R>::make_relation(std::string(concrete_relation), std::string(gov_view));
    }

    virtual const G make_gov_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      // look it up in the mapping dictionary
      DEBUG << "setting " << token.original() << " to " << word_mapper_.representation(token.original());
      return word_mapper_.representation(token.original());
    }

    virtual const G make_head_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      // look it up in the mapping dictionary
      return word_mapper_.representation(token.original());
    }
  };

  template <typename G, typename R>
  class SituationSuggestedPruner : public EntityMentionPruner<Mention<G, R> > {
  protected:
    typedef Vocabulary<G> GVocab;
    typedef Vocabulary<R> RVocab;
  private:
    // Record entity mention IDs to tokenizations
    const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization;
    // Record tokenizations to Situation mentions
    const concrete_util::uuid_map< std::list< concrete::SituationMention > > tok_id_to_sm;
    concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention;
     
    std::string dep_parse_toolname = "col-ccproc-deps";
    GVocab* gvocab_;
    RVocab* rvocab_;
    Vocabulary<std::string>* g_lemma_vocab_;
    Vocabulary<std::string>* r_lemma_vocab_;
    typedef Mention<G,R> M;
  protected:
    bool use_lexical_;
    bool add_to_vocabs_;
    StopWordList gov_stopwords_;
    bool use_gov_stopwords_;
  public:
    SituationSuggestedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc, 
				    Vocabulary<std::string>* g_lemma_vocab,
				    Vocabulary<std::string>* r_lemma_vocab,
				    bool use_lex) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(g_lemma_vocab), r_lemma_vocab_(r_lemma_vocab),
      use_lexical_(use_lex), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }
    SituationSuggestedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(false), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }
    SituationSuggestedPruner< G, R >(const concrete::Communication& comm,
				    GVocab *gvoc, RVocab* rvoc, bool use_lex) : 
    EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      gvocab_(gvoc), rvocab_(rvoc), g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(use_lex), add_to_vocabs_(true), use_gov_stopwords_(false) {
    }

    SituationSuggestedPruner< G, R >(const concrete::Communication& comm) : EntityMentionPruner< Mention< G, R> >(comm), 
      mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, "Stanford")),
      tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, "Semafor")),
      mention_id_to_mention(concrete_util::mention_id_to_mention(comm, "Semafor")), 
      g_lemma_vocab_(NULL), r_lemma_vocab_(NULL),
      use_lexical_(false), add_to_vocabs_(false), use_gov_stopwords_(false) {
    }

    void use_gov_stopwords(bool b) {
      use_gov_stopwords_ = b;
    }
    void gov_stopwords(const StopWordList& list_) {
      gov_stopwords_ = list_;
    }

    inline void gov_lemma_vocab(Vocabulary<std::string>* gv) {
      g_lemma_vocab_ = gv;
    }
    inline void rel_lemma_vocab(Vocabulary<std::string>* rv) {
      r_lemma_vocab_ = rv;
    }

    inline void use_lexical(bool b) {
      use_lexical_ = b;
    }
    bool use_lexical() {
      return use_lexical_;
    }
    bool use_lexical() const  {
      return use_lexical_;
    }
    virtual const R make_relation(const std::string rel, const std::string gov) const {
      std::string rel_ret = rel;
      rel_ret += "-" + gov;
      return rel_ret;
    }
    virtual const R make_relation(const std::string& concrete_relation, const std::string& dependency_rel,
				  G gov_view, const std::string& gov_lemma) const {
      return use_lexical_ ? make_relation(std::string(dependency_rel), std::string(gov_lemma)) : 
	make_relation(std::string(concrete_relation), std::string(gov_view));
    }
    virtual const G make_gov_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      return std::string(use_lexical_ ? token.lemma() : view_str);
    }
    virtual const G make_head_view(const AnnotatedToken< G >& token, const std::string& view_str) const {
      return std::string(use_lexical_ ? token.lemma() : view_str);
    }
    virtual void prune_without_sm(const concrete::EntityMention& conc_mention, 
				  const concrete::Tokenization& mention_tokenization,
				  int anchor_index,
				  std::vector<M>* created_mentions) const {
      //get the dependency parse
      const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(mention_tokenization, dep_parse_toolname);
      if(dep_parse == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " does not have any dependency parses with name containing " << dep_parse_toolname;
	return ;
      }
      //get the SituationMention corresponding to this mention's tokenization
      const concrete::TokenTagging * const posTagging = concrete_util::first_pos_tagging(mention_tokenization, "Stanford");
      if(posTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	return ;
      }
      const concrete::TokenTagging * const lemmaTagging = concrete_util::first_lemma_tagging(mention_tokenization, "Stanford");
      if(lemmaTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no LEMMA tagging";
	return ;
      }
      // for every dependency, check to see if dependency.dep == anchor_index
      for(const concrete::Dependency dependency : dep_parse->dependencyList) {
	if(dependency.dep != anchor_index || dependency.gov < 0) continue;
	// TODO: this really should be "generalizing"
	if(posTagging->taggedTokenList[dependency.gov].tag[0] != 'V') continue;
	if(use_gov_stopwords_ && 
	   gov_stopwords_.contains(lemmaTagging->taggedTokenList[dependency.gov].tag)) {
	  continue;
	}
	M mention;
	mention.latent(true);
	const std::string g_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	AnnotatedToken< G > gov(g_lemma,
				mention_tokenization.tokenList.tokenList[dependency.gov].text,
				posTagging->taggedTokenList[dependency.gov].tag);
	if(add_to_vocabs_) {
	  if(g_lemma_vocab_ != NULL) {
	    g_lemma_vocab_->make_word(g_lemma);
	  }
	}
	const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	AnnotatedToken< G > head(head_lemma,
				 mention_tokenization.tokenList.tokenList[dependency.dep].text,
				 posTagging->taggedTokenList[dependency.dep].tag);
	mention.gov(gov);
	mention.head(head);
	const std::string r_lemma = this->make_relation(dependency.edgeType, g_lemma);
	mention.rel_str( r_lemma );
	if(add_to_vocabs_) {
	  if(r_lemma_vocab_ != NULL) {
	    r_lemma_vocab_->make_word(r_lemma);
	  }
	}
	created_mentions->push_back(mention);
      }
    }
    virtual std::vector<M> prune(const concrete::EntityMention& conc_mention) const {
      const int anchor_index = conc_mention.tokens.anchorTokenIndex;
      std::vector<M> created_mentions;
      const concrete::Tokenization mention_tokenization = mention_id_to_tokenization.at(conc_mention.uuid);
      if(anchor_index < 0) {
	DEBUG << "Anchor index not set in Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
	// we should probably log a warning here
	return created_mentions;
      }
      if(tok_id_to_sm.find(mention_tokenization.uuid) == tok_id_to_sm.end()) {
	DEBUG << "No situation mentions found for Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
	prune_without_sm(conc_mention, mention_tokenization, anchor_index, &created_mentions);
	return created_mentions;
      }

      //get the dependency parse
      const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(mention_tokenization, dep_parse_toolname);
      if(dep_parse == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " does not have any dependency parses with name containing " << dep_parse_toolname;
	return created_mentions;
      }
      //get the SituationMention corresponding to this mention's tokenization
      std::list<concrete::SituationMention> sit_ment_list = tok_id_to_sm.at(mention_tokenization.uuid);
      const concrete::TokenTagging * const posTagging = concrete_util::first_pos_tagging(mention_tokenization, "Stanford");
      if(posTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	return created_mentions;
      }
      const concrete::TokenTagging * const lemmaTagging = concrete_util::first_lemma_tagging(mention_tokenization, "Stanford");
      if(lemmaTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no LEMMA tagging";
	return created_mentions;
      }
      // for every dependency, check to see if dependency.dep == anchor_index
      for(const concrete::Dependency dependency : dep_parse->dependencyList) {
	if(dependency.dep != anchor_index || dependency.gov < 0) continue;
	// TODO: this really should be "generalizing"
	if(posTagging->taggedTokenList[dependency.gov].tag[0] != 'V') continue;
	if(use_gov_stopwords_ && 
	   gov_stopwords_.contains(lemmaTagging->taggedTokenList[dependency.gov].tag)) {
	  continue;
	}
	// now verify that the governor corresponds to a situation mention trigger
	for(const concrete::SituationMention sit_ment : sit_ment_list) {
	  if(situation_mention_violations(sit_ment) || 
	     !concrete::util::index_overlap_mention<concrete::SituationMention>(sit_ment, dependency.gov)) continue;
	  // now, make sure that the dependency token is covered by one of the arcs
	  int mention_arg_idx = -1;
	  int mention_counter = -1;
	  for(const concrete::MentionArgument mention_argument : sit_ment.argumentList) {
	    ++mention_counter;
	    if(!mention_argument.__isset.role) {
	      WARN << "Skipping mention argument " << mention_counter << " due to missing `role` field";
	      continue;
	    }
	    if(!mention_argument.__isset.entityMentionId) {
	      WARN << "Skipping mention argument " << mention_counter << " due to missing `entityMentionId` field (this object only works with EntityMention arguments)";
	      continue;
	    }
	    // const concrete::Tokenization arg_mention_tokenization = 
	    //   sm_arg_mention_id_to_tokenization.at(mention_argument.entityMentionId);
	    const concrete::EntityMention arg_ent_mention =
	      mention_id_to_mention.at(mention_argument.entityMentionId);
	    if(!concrete::util::index_overlap_mention<concrete::EntityMention>(arg_ent_mention, dependency.dep)) continue;
	    mention_arg_idx = mention_counter;
	    break;
	  }
	  M mention;
	  if(mention_arg_idx < 0) {
	    mention.latent(true);
	    const std::string g_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	    AnnotatedToken< G > gov(g_lemma,
				    mention_tokenization.tokenList.tokenList[dependency.gov].text,
				    posTagging->taggedTokenList[dependency.gov].tag);
	    if(add_to_vocabs_) {
	      if(g_lemma_vocab_ != NULL) {
		g_lemma_vocab_->make_word(g_lemma);
	      }
	    }
	    const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	    AnnotatedToken< G > head(head_lemma,
				     mention_tokenization.tokenList.tokenList[dependency.dep].text,
				     posTagging->taggedTokenList[dependency.dep].tag);
	    mention.gov(gov);
	    mention.head(head);
	    const std::string r_lemma = this->make_relation(dependency.edgeType, g_lemma);
	    mention.rel_str( r_lemma );
	    if(add_to_vocabs_) {
	      if(r_lemma_vocab_ != NULL) {
		r_lemma_vocab_->make_word(r_lemma);
	      }
	    }
	  } else {
	    mention.latent(false);
	    const concrete::MentionArgument rel_mention_arg = sit_ment.argumentList[mention_arg_idx];
	    const std::string g_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	    AnnotatedToken< G > gov(g_lemma,
				    mention_tokenization.tokenList.tokenList[dependency.gov].text,
				    posTagging->taggedTokenList[dependency.gov].tag);
	    // The gov view is 
	    const std::string frame_trigger = sit_ment.situationKind;
	    const G gov_view = this->make_gov_view(gov, frame_trigger);
	    gov.view( gov_view );
	    if(add_to_vocabs_) {
	      gvocab_->make_word(gov_view);
	      if(g_lemma_vocab_ != NULL) {
		g_lemma_vocab_->make_word(g_lemma);
	      }
	    }
	    const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	    AnnotatedToken< G > head(head_lemma,
				     mention_tokenization.tokenList.tokenList[dependency.dep].text,
				     posTagging->taggedTokenList[dependency.dep].tag);
	    head.view( this->make_head_view(head, head_lemma) );
	    mention.gov(gov);
	    mention.head(head);
	    const R rel = this->make_relation(rel_mention_arg.role, dependency.edgeType, gov.view(), g_lemma);
	    const std::string r_lemma = this->make_relation(dependency.edgeType, g_lemma);
	    mention.rel_str( r_lemma );
	    mention.rel( rel );
	    if(add_to_vocabs_) {
	      rvocab_->make_word(rel);
	      if(r_lemma_vocab_ != NULL) {
		r_lemma_vocab_->make_word(r_lemma);
	      }
	    }
	  }
	  created_mentions.push_back(mention);
	}
      }
      return created_mentions;
    }

    /**
     * Return the number of assumption violations.
     * Any return value > 0 means we will skip processing this mention.
     */
    int situation_mention_violations(const concrete::SituationMention& sm) const {
      int num = 0;
      if(! sm.__isset.tokens) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to missing 'tokens' field";
	++num;
      }
      if(! sm.__isset.situationKind) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to missing 'situationKind' field";
	++num;
      }
      if(!sm.argumentList.size()) {
	DEBUG << "Skipping situation mention " << sm.uuid.uuidString << " due to empty 'mentionArgument' field";
	++num;
      }
      return num;
    }
  };

}

#endif
