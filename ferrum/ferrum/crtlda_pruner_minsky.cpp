#include "ferrum/concrete.hpp"
#include "ferrum/crtlda_pruner_minsky.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/tm_pruner.hpp"

namespace ferrum {
  std::istream& operator>>(std::istream& in, ferrum::MinskyPrunerEnum &how) {
    std::string token;
    in >> token;
    if (token == "VerbGoverned")
      how = ferrum::MinskyPrunerEnum::VerbGoverned;
    else if (token == "SituationGoverned" )
      how = ferrum::MinskyPrunerEnum::SituationGoverned;
    else if (token == "DocBOW")
      how = ferrum::MinskyPrunerEnum::DocBOW;
    else if (token == "VerbBOW")
      how = ferrum::MinskyPrunerEnum::VerbBOW;
    else {
      ERROR << "Unknown pruner choice " << token;
      throw 5;
    }
    //else throw boost::program_options::validation_error("Invalid unit");
    return in;
  }
  std::ostream& operator<<(std::ostream& out, ferrum::MinskyPrunerEnum &how) {
    switch(how) {
    case ferrum::MinskyPrunerEnum::VerbGoverned:
      out << "VerbGoverned";
      break;
    case ferrum::MinskyPrunerEnum::SituationGoverned:
      out << "SituationGoverned";
      break;
    case ferrum::MinskyPrunerEnum::DocBOW:
      out << "DocBOW";
    case ferrum::MinskyPrunerEnum::VerbBOW:
      out << "VerbBOW";
      break;
    default:
      ERROR << "Unknown pruner choice " << how;
      throw 5;
    }
    //else throw boost::program_options::validation_error("Invalid unit");
    return out;
  }
  MinskyVerbGovernedPruner::MinskyVerbGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* g_lemma_vocab,
   Vocabulary<std::string>* r_lemma_vocab,
   const ferrum::Toolnames& tools
   ) :
    MinskyPruner<MinskyVerbGovernedPruner>(comm, tools),
    mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, tools.entity_mention_tool)),
    g_lemma_vocab_(g_lemma_vocab),
    r_lemma_vocab_(r_lemma_vocab),
    ndh_(1) {
    dep_parse_toolname = tools.dep_parse_tool;
  }
  MinskyVerbGovernedPruner::MinskyVerbGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* g_lemma_vocab,
   Vocabulary<std::string>* r_lemma_vocab
   ) :
    MinskyVerbGovernedPruner(comm, g_lemma_vocab, r_lemma_vocab, Toolnames()) {
  }
  void MinskyVerbGovernedPruner::num_dep_hops(unsigned int ndh) {
    ndh_ = ndh;
  }
  std::string MinskyVerbGovernedPruner::make_relation(const std::string& concrete_relation, const std::string& gov_view) const {
    std::string ret_str(concrete_relation + "-" + gov_view);
    return ret_str;
  }
  std::string MinskyVerbGovernedPruner::make_gov_view(const std::string& lemma) const {
    return lemma;
  }
  std::string MinskyVerbGovernedPruner::make_head_view(const std::string& lemma) const {
    return lemma;
  }
  const concrete::Tokenization&
  MinskyVerbGovernedPruner::mention_tokenization(const concrete::UUID& uuid) {
    return mention_id_to_tokenization.at(uuid);
  }
  std::vector<minsky::Mention> MinskyVerbGovernedPruner::prune(const concrete::EntityMention& conc_mention) const {
    const concrete::Tokenization& mention_tokenization = mention_id_to_tokenization.at(conc_mention.uuid);
    const int anchor_index = conc_mention.tokens.anchorTokenIndex;
    std::vector<M> created_mentions;
    if(anchor_index < 0) {
      ERROR << "Anchor index not set in Concrete entity mention " << conc_mention.uuid.uuidString << ", in Tokenization " << mention_tokenization.uuid.uuidString;
      return created_mentions;
    }
    //get the dependency parse
    const concrete::DependencyParse* const dep_parse = concrete_util::first_dependency_parse(mention_tokenization, dep_parse_toolname);
    if(dep_parse == NULL) {
      WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " does not have any dependency parses with name containing " << dep_parse_toolname;
      return created_mentions;
    }
    // for every dependency, check to see if dependency.dep == anchor_index
    int dep_id = -1;
    for(const concrete::Dependency& dependency : dep_parse->dependencyList) {
      ++dep_id;
      if(dependency.dep != anchor_index || dependency.gov < 0) continue;
      const concrete::TokenTagging * const posTagging = concrete_util::first_pos_tagging(mention_tokenization, "Stanford");
      if(posTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	continue;
      }
      PathGeneralizer pg(&(posTagging->taggedTokenList), &(dep_parse->dependencyList), ndh_);
      if(! pg.generalize_to_v( (size_t)dep_id, 0) ) {
	continue;
      }
      //if(posTagging->taggedTokenList[dependency.gov].tag[0] != 'V') continue;
      const concrete::TokenTagging * const lemmaTagging = concrete_util::first_lemma_tagging(mention_tokenization, "Stanford");
      if(lemmaTagging == NULL) {
	WARN << "Tokenization " << mention_tokenization.uuid.uuidString << " has no POS tagging";
	continue;
      }
      M mention;
      mention.__set_location_id(mention_tokenization.uuid.uuidString);
      //mention.latent(false);

      minsky::PredArg syntactic_pa;
      syntactic_pa.__set_annot_level(minsky::AnnotationLevel::SYNTAX);
      minsky::RelationFiller pred;
      //const std::string& gov_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
      std::string gov_lemma = pg.form_gov_tagging(lemmaTagging->taggedTokenList);
      ferrum::lower(gov_lemma);
      pred.__set_word((*g_lemma_vocab_)(gov_lemma));
      pred.__set_position(dependency.gov);
      syntactic_pa.__set_predicate(pred);
      std::string r_lemma1 = pg.form_arc_tagging();
      std::string r_lemma = this->make_relation(r_lemma1, gov_lemma);
      ferrum::lower(r_lemma);
      //std::string r_lemma = this->make_relation(dependency.edgeType, gov_lemma);
      syntactic_pa.__set_relation((*r_lemma_vocab_)(r_lemma));
      //INFO << r_lemma << "(" << gov_lemma << ", " << pg.form_dep_tagging(lemmaTagging->taggedTokenList) << ")";
      // if(pg.size() >= 3) {
      // 	const std::list<GeneralizedDependencyPath>& al = pg.get_path();
      // 	for(const GeneralizedDependencyPath& gdp: al) {
      // 	  std::cout << gdp.dependency_index << " ";
      // 	}
      // 	std::cout << std::endl;
      // }
      mention.structures.push_back(syntactic_pa);
      mention.__isset.structures = true;
      mention.__set_id(conc_mention.uuid.uuidString);
      created_mentions.push_back(mention);
    }
    return created_mentions;
  }

  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string> *gvoc, Vocabulary<std::string>* rvoc, 
   Vocabulary<std::string>* g_lemma_vocab,
   Vocabulary<std::string>* r_lemma_vocab,
   const ferrum::Toolnames& tools,
   bool use_lex) :
    MinskyPruner<MinskySituationGovernedPruner>(comm, tools),
    mention_id_to_tokenization(concrete_util::mention_id_to_tokenization(comm, tools.entity_mention_tool)),
    tok_id_to_sm(concrete_util::tokenization_id_to_situation_mention(comm, tools.situation_mention_tool)),
    mention_id_to_mention(concrete_util::mention_id_to_mention(comm, tools.situation_mention_tool)), 
    dep_parse_toolname(tools.dep_parse_tool),
    gvocab_(gvoc),
    rvocab_(rvoc),
    g_lemma_vocab_(g_lemma_vocab),
    r_lemma_vocab_(r_lemma_vocab),
    add_syntax_(false),
    add_semantics_(false),
    use_lexical_(use_lex),
    add_to_vocabs_(true),
    use_gov_stopwords_(false) {
    if(gvocab_ && rvocab_) {
      add_semantics_ = true;
    }
    if(g_lemma_vocab_ && r_lemma_vocab_) {
      add_syntax_ = true;
    }
  }
  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string> *gvoc, Vocabulary<std::string>* rvoc, 
   Vocabulary<std::string>* g_lemma_vocab,
   Vocabulary<std::string>* r_lemma_vocab,
   bool use_lex) : 
    MinskySituationGovernedPruner(comm,
				  gvoc,
				  rvoc,
				  g_lemma_vocab,
				  r_lemma_vocab,
				  ferrum::Toolnames(),
				  false) {
  }
  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string> *gvoc, Vocabulary<std::string>* rvoc) :
    MinskySituationGovernedPruner(comm,
				  NULL,
				  NULL,
				  gvoc,
				  rvoc,
				  false) {
  }
  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string> *gvoc, Vocabulary<std::string>* rvoc,
   const ferrum::Toolnames& tools) :
    MinskySituationGovernedPruner(comm,
				  NULL,
				  NULL,
				  gvoc,
				  rvoc,
				  tools,
				  false) {
  }
  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string> *gvoc,
   Vocabulary<std::string>* rvoc,
   bool use_lex
   ) : 
    MinskySituationGovernedPruner(comm, NULL, NULL, gvoc, rvoc, use_lex) {
  }
  MinskySituationGovernedPruner::MinskySituationGovernedPruner
  (
   const concrete::Communication& comm
   ) :
    MinskySituationGovernedPruner(comm, NULL, NULL, NULL, NULL, false) {
  }
  std::string MinskySituationGovernedPruner::make_relation(const std::string rel, const std::string gov) const {
    std::string rel_ret = rel;
    rel_ret += "-" + gov;
    return rel_ret;
  }
  std::string MinskySituationGovernedPruner::make_relation(const std::string& concrete_relation, const std::string& dependency_rel,
							   const std::string& gov_view, const std::string& gov_lemma) const {
    return use_lexical_ ?
      make_relation(dependency_rel, gov_lemma) :
      make_relation(concrete_relation, gov_view);
  }
  std::string MinskySituationGovernedPruner::make_gov_view(const std::string& lemma, const std::string& view_str) const {
    return use_lexical_ ? lemma : view_str;
  }
  std::string MinskySituationGovernedPruner::make_head_view(const std::string& lemma, const std::string& view_str) const {
    return use_lexical_ ? lemma : view_str;
  }

  const concrete::Tokenization&
  MinskySituationGovernedPruner::mention_tokenization(const concrete::UUID& uuid) {
    return mention_id_to_tokenization.at(uuid);
  }

  std::vector<minsky::Mention> MinskySituationGovernedPruner::prune
  (
   const concrete::EntityMention& conc_mention
   ) const {
    const concrete::Tokenization mention_tokenization = mention_id_to_tokenization.at(conc_mention.uuid);
    const int anchor_index = conc_mention.tokens.anchorTokenIndex;
    std::vector<minsky::Mention> created_mentions;
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
    for(const concrete::Dependency& dependency : dep_parse->dependencyList) {
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
	   !(concrete::util::index_overlap_mention<concrete::SituationMention>(sit_ment, dependency.gov))) continue;
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
	minsky::Mention mention;
	std::vector<minsky::PredArg> mention_predargs;
	std::string g_lemma = lemmaTagging->taggedTokenList[dependency.gov].tag;
	ferrum::lower(g_lemma);
	if(add_syntax_) {
	  minsky::PredArg syntactic_pa;
	  syntactic_pa.__set_annot_level(minsky::AnnotationLevel::SYNTAX);
	  minsky::RelationFiller pred;
	  pred.__set_word((*g_lemma_vocab_)(g_lemma));
	  pred.__set_position(dependency.gov);
	  syntactic_pa.__set_predicate(pred);
	  std::string r_lemma = this->make_relation(dependency.edgeType, g_lemma);
	  ferrum::lower(r_lemma);
	  syntactic_pa.__set_relation((*r_lemma_vocab_)(r_lemma));
	  mention_predargs.push_back(syntactic_pa);
	}
	if(add_semantics_) {
	  minsky::PredArg semantic_pa;
	  semantic_pa.__set_annot_level(minsky::AnnotationLevel::SEMANTIC);
	  minsky::RelationFiller pred;
	  const std::string frame_trigger = sit_ment.situationKind;
	  std::string gov_view = this->make_gov_view(g_lemma, frame_trigger);
	  pred.__set_word((*gvocab_)(gov_view));
	  pred.__set_position(dependency.gov);
	  semantic_pa.__set_predicate(pred);
	  std::string rel = this->make_relation(rel_mention_arg.role, dependency.edgeType, gov_view, g_lemma);
	  semantic_pa.__set_relation((*rvocab_)(rel));
	  mention_predargs.push_back(semantic_pa);
	}
	// const std::string head_lemma = lemmaTagging->taggedTokenList[dependency.dep].tag;
	// AnnotatedToken< G > head(head_lemma,
	// 			   mention_tokenization.tokenList.tokenList[dependency.dep].text,
	// 			   posTagging->taggedTokenList[dependency.dep].tag);
	// head.view( this->make_head_view(head, head_lemma) );
	if(mention_predargs.size()) {
	  mention.__set_structures(mention_predargs);
	}
	mention.__set_id(conc_mention.uuid.uuidString);
	created_mentions.push_back(mention);
      } // end loop over Situation Mentions
    } // end loop over each dependency relation
    return created_mentions;
  }

    
  /**
   * Return the number of assumption violations.
   * Any return value > 0 means we will skip processing this mention.
   */
  int MinskySituationGovernedPruner::situation_mention_violations(const concrete::SituationMention& sm) const {
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

  ////////////////////////////////////////////////

  MinskyDocBOWPruner::MinskyDocBOWPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* vocab,
   minsky::WordAnnotation::type wannot,
   unsigned int
   ) :
    MinskyDocBOWPruner(comm, vocab, ferrum::Toolnames(), wannot) {
  }
  MinskyDocBOWPruner::MinskyDocBOWPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* vocab,
   const ferrum::Toolnames& tools,
   minsky::WordAnnotation::type wannot,
   unsigned int
   ) :
    MinskySentencePruner<MinskyDocBOWPruner, concrete::Communication>(comm, tools, wannot),
    vocab_(vocab) {
    dep_parse_toolname = tools.dep_parse_tool;
  }

  minsky::WordsClause
  MinskyDocBOWPruner::clause_create(const concrete::Communication& comm) const {
    minsky::WordsClause wc;
    minsky::CountList cl;
    switch(annot_) {
    case minsky::WordAnnotation::ORTHOGRAPHIC:
      for(const concrete::Token& token : concrete::util::TokenIterator(comm)) {
	if(! token.__isset.text) {
	  ERROR << "Token.text needs to be set";
	  continue;
	}
	std::string ltext = token.text;
	ferrum::lower(ltext);
	cl.icounts[ vocab_->operator()(ltext) ]++;
      }
      break;
    case minsky::WordAnnotation::LEMMA:
      for(const concrete::Tokenization& tkz : concrete::util::TokenizationIterator(comm)) {
	const concrete::TokenTagging* const lemmas = concrete_util::first_set(tkz.tokenTaggingList, "LEMMA");
	if(lemmas == NULL) continue;
	for(const concrete::TaggedToken& tt : lemmas->taggedTokenList) {
	  if(! tt.__isset.tag) {
	    ERROR << "TaggedToken.tag needs to be set";
	    continue;
	  }
	  std::string ltag = tt.tag;
	  ferrum::lower(ltag);
	  cl.icounts[ vocab_->operator()(tt.tag) ]++;
	}
      }
      break;
    default:
      ERROR << "Unknown type " << annot_;
      break;
    }
    if(cl.icounts.size()) {
      cl.__isset.icounts = true;
      wc.__set_counts(cl);
    }
    return wc;
  }

  ////////////////////////////////////////////////

  MinskyVerbBOWPruner::MinskyVerbBOWPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* vocab,
   minsky::WordAnnotation::type wannot,
   unsigned int ndh
   ) :
    MinskyVerbBOWPruner(comm, vocab, ferrum::Toolnames(), wannot, ndh) {
  }
  MinskyVerbBOWPruner::MinskyVerbBOWPruner
  (
   const concrete::Communication& comm,
   Vocabulary<std::string>* vocab,
   const ferrum::Toolnames& tools,
   minsky::WordAnnotation::type wannot,
   unsigned int ndh
   ) :
    MinskySentencePruner<MinskyVerbBOWPruner, concrete::Communication>(comm, tools, wannot),
    vocab_(vocab),
    tools_(tools),
    ndh_(ndh) {
    dep_parse_toolname = tools_.dep_parse_tool;
  }

  minsky::WordsClause
  MinskyVerbBOWPruner::clause_create(const concrete::Communication& comm) const {
    minsky::WordsClause wc;
    minsky::CountList cl;
    switch(annot_) {
    case minsky::WordAnnotation::ORTHOGRAPHIC:
      // for(const concrete::Token& token : concrete::util::TokenIterator(comm)) {
      // 	if(! token.__isset.text) {
      // 	  ERROR << "Token.text needs to be set";
      // 	  continue;
      // 	}
      // 	cl.icounts[ vocab_->operator()(token.text) ]++;
      // }
      //break;
    case minsky::WordAnnotation::LEMMA:
      {
	ferrum::Vocabulary<std::string> rel_vocab;
	minsky::EDoc edoc = ferrum::MinskyVerbGovernedPruner::make(comm, ndh_, vocab_, &rel_vocab, tools_);
	for(const minsky::Entity& entity : edoc.entities) {
	  for(const minsky::Mention& mention : entity.mentions) {
	    for(const minsky::PredArg& pa : mention.structures) {
	      if(pa.annot_level != minsky::AnnotationLevel::SYNTAX) continue;
	      if(pa.__isset.predicate && pa.predicate.__isset.word) {
		cl.icounts[ pa.predicate.word ]++;
	      }
	    }
	  }
	}
      }
      break;
    default:
      ERROR << "Unknown type " << annot_;
      break;
    }
    if(cl.icounts.size()) {
      cl.__isset.icounts = true;
      wc.__set_counts(cl);
    }
    return wc;
  }
}


