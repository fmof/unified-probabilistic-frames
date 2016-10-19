#ifndef FERRUM_SENTENCES_HPP_
#define FERRUM_SENTENCES_HPP_

#include "ferrum/concrete.hpp"

namespace ferrum {
  
  class FlatTokenization {
  public:
    enum Form {
      UNSPECIFIED = -1,
      ORTHO = 0,
      LEMMA = 1,
      POS = 2,
      NER = 3
    };
    FlatTokenization(const concrete::Tokenization& tokenization, Form form = UNSPECIFIED);
    const std::string& operator()(size_t i) const;
    const std::string& word(size_t i) const;
    const std::string& lemma(size_t i) const;
    const std::string& pos(size_t i) const;
    const std::string& ner(size_t i) const;

    bool empty_lemma() const;
    bool empty_pos() const;
    bool empty_ner() const;

    bool is_verb(size_t i) const;
    void form(Form f);
  private:
    std::vector<concrete::Token> words_;
    std::vector<concrete::TaggedToken> lemma_;
    std::vector<concrete::TaggedToken> pos_;
    std::vector<concrete::TaggedToken> ner_;
    Form form_;
  };

  class DependencyGraph {
  public:
    DependencyGraph(size_t num_words);
    const std::string& edge(int, int) const;

    /**
     * Predicate is a function (or functor) of type
     *   bool operator()(size_t) const
     * It must operate on 
     */
    template <typename Predicate>
    DependencyGraph& fill(const concrete::DependencyParse& dp, Predicate pred, bool inverse);
    size_t num_roots() const;
    const std::set<int>& roots() const;

    template <typename ReturnOp, typename NextOp>
    void bfs_from(int, ReturnOp&, const NextOp&) const;
    template <typename ReturnOp>
    void bfs_from(int, ReturnOp&) const;

    template <typename ReturnOp, typename NextOp>
    void bfs(ReturnOp&, const NextOp&) const;
  private:
    size_t num_vertices_;
    std::set<int> roots_;
    std::vector< std::list<int> > graph_;
    std::vector< std::vector<std::string> > edges_;
  };

  class WordChain {
  public:
    WordChain(const FlatTokenization&);
    void operator()(size_t);
    bool operator()(int, int, const DependencyGraph&) const;
    const std::vector<std::string>& get() const;
  private:
    FlatTokenization ft_;
    std::vector<std::string> chain_;
  };
}

#include "ferrum/sentences.tcc"

#endif
