#ifndef FERRUM_LIBNAR_CRTLDA_MINSKY_PRUNER_H_
#define FERRUM_LIBNAR_CRTLDA_MINSKY_PRUNER_H_

#include "ferrum/util.hpp"
#include "ferrum/crtlda_defs.hpp"
#include "ferrum/crtlda_pruner_base.hpp"
#include "ferrum/tm_pruner.hpp"
#include "ferrum/concrete.hpp"
#include "ferrum/minsky.hpp"
#include "ferrum/words.hpp"

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
  class MinskyPruner : public EntityMentionPruner< minsky::Mention > {
  protected:
    ferrum::Toolnames tools_;
    void make_edoc(minsky::EDoc* doc);
    virtual const concrete::Tokenization& mention_tokenization(const concrete::UUID&) = 0;
  public:
    MinskyPruner<Subpruner>(const concrete::Communication& comm,
			    const ferrum::Toolnames tools);
    virtual ~MinskyPruner();
    template <typename ...ArgTypes>
    static minsky::EDoc make(const concrete::Communication& comm,
			     unsigned int num_dep_hops,
			     ArgTypes... args);
    template <typename ...ArgTypes>
    static minsky::EDoc make(const concrete::Communication& comm,
			     ArgTypes... args);
    virtual void num_dep_hops(unsigned int ndh);
  };

  class MinskyVerbGovernedPruner : public MinskyPruner< MinskyVerbGovernedPruner > {
  private:
    const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization;
    std::string dep_parse_toolname = "col-ccproc-deps";
    Vocabulary<std::string>* g_lemma_vocab_;
    Vocabulary<std::string>* r_lemma_vocab_;
    unsigned int ndh_;
    typedef minsky::Mention M;
  protected:
    virtual const concrete::Tokenization& mention_tokenization(const concrete::UUID&);
  public:
    MinskyVerbGovernedPruner(const concrete::Communication& comm, Vocabulary<std::string>* g_lemma_vocab, Vocabulary<std::string>* r_lemma_vocab);
    MinskyVerbGovernedPruner(const concrete::Communication& comm, Vocabulary<std::string>* g_lemma_vocab, Vocabulary<std::string>* r_lemma_vocab, const Toolnames& tools);
    virtual std::string make_relation(const std::string& concrete_relation, const std::string& gov_view) const;
    virtual std::string make_gov_view(const std::string& lemma) const;
    virtual std::string make_head_view(const std::string& lemma) const;
    virtual std::vector<minsky::Mention> prune(const concrete::EntityMention& conc_mention) const;
    virtual void num_dep_hops(unsigned int ndh);
  };

  class MinskySituationGovernedPruner : public MinskyPruner<MinskySituationGovernedPruner> {
  private:
    // Record entity mention IDs to tokenizations
    const concrete_util::uuid_map<concrete::Tokenization> mention_id_to_tokenization;
    // Record tokenizations to Situation mentions
    const concrete_util::uuid_map< std::list< concrete::SituationMention > > tok_id_to_sm;
    concrete_util::uuid_map<concrete::EntityMention> mention_id_to_mention;
     
    std::string dep_parse_toolname;
    Vocabulary<std::string>* gvocab_;
    Vocabulary<std::string>* rvocab_;
    Vocabulary<std::string>* g_lemma_vocab_;
    Vocabulary<std::string>* r_lemma_vocab_;
    bool add_syntax_;
    bool add_semantics_;
  protected:
    bool use_lexical_;
    bool add_to_vocabs_;
    StopWordList gov_stopwords_;
    bool use_gov_stopwords_;
    virtual const concrete::Tokenization& mention_tokenization(const concrete::UUID&);
  public:
    MinskySituationGovernedPruner(const concrete::Communication& comm,
				  Vocabulary<std::string> *gvoc,
				  Vocabulary<std::string>* rvoc,
				  Vocabulary<std::string>* g_lemma_vocab,
				  Vocabulary<std::string>* r_lemma_vocab,
				  const ferrum::Toolnames& tools,
				  bool use_lex);
    MinskySituationGovernedPruner(const concrete::Communication& comm,
				  Vocabulary<std::string> *gvoc,
				  Vocabulary<std::string>* rvoc,
				  Vocabulary<std::string>* g_lemma_vocab,
				  Vocabulary<std::string>* r_lemma_vocab,
				  bool use_lex);
    MinskySituationGovernedPruner(const concrete::Communication& comm,
				  Vocabulary<std::string> *gvoc,
				  Vocabulary<std::string>* rvoc);
    MinskySituationGovernedPruner(const concrete::Communication& comm,
				  Vocabulary<std::string> *gvoc,
				  Vocabulary<std::string>* rvoc,
				  const ferrum::Toolnames& tools);
    MinskySituationGovernedPruner(const concrete::Communication& comm,
				  Vocabulary<std::string> *gvoc,
				  Vocabulary<std::string>* rvoc,
				  bool use_lex);
    MinskySituationGovernedPruner(const concrete::Communication& comm);

    inline void use_gov_stopwords(bool b) {
      use_gov_stopwords_ = b;
    }
    inline void gov_stopwords(const StopWordList& list_) {
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
    inline bool use_lexical() {
      return use_lexical_;
    }
    inline bool use_lexical() const  {
      return use_lexical_;
    }
    virtual std::string make_relation(const std::string rel, const std::string gov) const;
    virtual std::string make_relation(const std::string& concrete_relation, const std::string& dependency_rel,
				      const std::string& gov_view, const std::string& gov_lemma) const;
    virtual std::string make_gov_view(const std::string& token, const std::string& view_str) const;
    virtual std::string make_head_view(const std::string& token, const std::string& view_str) const;
    virtual std::vector<minsky::Mention> prune(const concrete::EntityMention& conc_mention) const;
    
    /**
     * Return the number of assumption violations.
     * Any return value > 0 means we will skip processing this mention.
     */
    int situation_mention_violations(const concrete::SituationMention& sm) const;
  };

  /////////////////////////////////////////////////////////////////////////

  template <typename Subpruner, typename ConcreteType>
  class MinskySentencePruner : public BOWPruner< minsky::WordsClause, ConcreteType > {
  protected:
    ferrum::Toolnames tools_;
    void make_simple_doc(minsky::SimpleDoc* doc);
    minsky::WordAnnotation::type annot_;
    
  public:
    MinskySentencePruner<Subpruner, ConcreteType>(const concrete::Communication& comm,
						  const ferrum::Toolnames tools,
						  minsky::WordAnnotation::type wannot);
    virtual ~MinskySentencePruner();
    template <typename ...ArgTypes>
    static minsky::SimpleDoc make(const concrete::Communication& comm,
				  unsigned int num_dep_hops,
				  ArgTypes... args);
    template <typename ...ArgTypes>
    static minsky::SimpleDoc make(const concrete::Communication& comm,
				  ArgTypes... args);
  };

  class MinskyDocBOWPruner : public MinskySentencePruner< MinskyDocBOWPruner, concrete::Communication > {
  private:
    std::string dep_parse_toolname = "col-ccproc-deps";
    Vocabulary<std::string>* vocab_;
  public:
    MinskyDocBOWPruner(const concrete::Communication& comm, Vocabulary<std::string>* vocab, minsky::WordAnnotation::type wannot = minsky::WordAnnotation::ORTHOGRAPHIC, unsigned int x = 2);
    MinskyDocBOWPruner(const concrete::Communication& comm, Vocabulary<std::string>* vocab, const Toolnames& tools, minsky::WordAnnotation::type wannot = minsky::WordAnnotation::ORTHOGRAPHIC, unsigned int x = 2);
    virtual minsky::WordsClause clause_create(const concrete::Communication& conc_obj) const;
  };

  ////////////////////////////////////////////////

  class MinskyVerbBOWPruner : public MinskySentencePruner< MinskyVerbBOWPruner, concrete::Communication > {
  private:
    std::string dep_parse_toolname = "col-ccproc-deps";
    Vocabulary<std::string>* vocab_;
    Toolnames tools_;
    unsigned int ndh_;
  public:
    MinskyVerbBOWPruner(const concrete::Communication& comm, Vocabulary<std::string>* vocab, minsky::WordAnnotation::type wannot = minsky::WordAnnotation::ORTHOGRAPHIC, unsigned int x = 2);
    MinskyVerbBOWPruner(const concrete::Communication& comm, Vocabulary<std::string>* vocab, const Toolnames& tools, minsky::WordAnnotation::type wannot = minsky::WordAnnotation::ORTHOGRAPHIC, unsigned int x = 2);
    virtual minsky::WordsClause clause_create(const concrete::Communication& conc_obj) const;
  };

  /////////////////////////////////////////////////////////////////

  enum MinskyPrunerEnum {
    UNSET = 0,
    VerbGoverned = 1,
    SituationGoverned = 2,
    DocBOW = 10,
    VerbBOW = 11
  };

  std::istream& operator>>(std::istream& in, ferrum::MinskyPrunerEnum &how);
  std::ostream& operator<<(std::ostream& in, ferrum::MinskyPrunerEnum &how);
} // end namespace crtlda

#include "ferrum/crtlda_pruner_minsky.tcc"

#endif
