#include "ferrum/sentences.hpp"

#include <boost/algorithm/string/case_conv.hpp>

namespace ferrum {
  
  FlatTokenization::FlatTokenization(const concrete::Tokenization& tokenization, Form f) {
    this->form(f);
    words_ = tokenization.tokenList.tokenList;
    typedef std::vector<concrete::TaggedToken>* ttptr;
    const size_t num = 3;
    ttptr fillees[num] = {&lemma_, &pos_, &ner_};
    std::string types[num] = {"LEMMA", "POS", "NER"};
    for(size_t i = 0; i < num; ++i) {
      const concrete::TokenTagging* const ptr = concrete_util::first_set(tokenization.tokenTaggingList, types[i]);
      if(ptr != NULL) {
	*(fillees[i]) = ptr->taggedTokenList;
      }
    }
  }
  const std::string& FlatTokenization::operator()(size_t i) const {
    switch(form_) {
    case ORTHO:
      return word(i);
    case LEMMA:
      return lemma(i);
    case POS:
      return pos(i);
    case NER:
      return ner(i);
    case UNSPECIFIED:
    default:
      ERROR << "Unknown form " << form_;
      throw 5;
    }
  }
  void FlatTokenization::form(Form f) {
    switch(f) {
    case ORTHO:
    case LEMMA:
    case POS:
    case NER:
      form_ = f;
      break;
    case UNSPECIFIED:
    default:
      ERROR << "Unknown form " << form_;
      throw 5;
    }
  }
  const std::string& FlatTokenization::word(size_t i) const {
    return words_.at(i).text;
  }
  const std::string& FlatTokenization::lemma(size_t i) const {
    return lemma_.at(i).tag;
  }
  const std::string& FlatTokenization::pos(size_t i) const {
    return pos_.at(i).tag;
  }
  const std::string& FlatTokenization::ner(size_t i) const {
    return ner_.at(i).tag;
  }
  bool FlatTokenization::empty_lemma() const {
    return lemma_.empty();
  }
  bool FlatTokenization::empty_pos() const {
    return pos_.empty();
  }
  bool FlatTokenization::empty_ner() const {
    return ner_.empty();
  }
  bool FlatTokenization::is_verb(size_t i) const {
    return pos(i)[0] == 'V';
  }

  DependencyGraph::DependencyGraph(size_t num_words) :
    num_vertices_(num_words + 1),
    roots_(),
    graph_(num_vertices_, std::list<int>()),
    edges_(num_vertices_, std::vector<std::string>(num_vertices_)) {     
  }
  size_t DependencyGraph::num_roots() const {
    return roots_.size();
  }
  const std::set<int>& DependencyGraph::roots() const {
    return roots_;
  }

  const std::string& DependencyGraph::edge(int p, int c) const {
    return edges_.at(p).at(c);
  }

  WordChain::WordChain(const FlatTokenization& ft) : ft_(ft) {
  }
  void WordChain::operator()(size_t i) {
    std::string word_str(ft_(i));
    boost::algorithm::to_lower(word_str);
    chain_.push_back(word_str);
  }
  bool WordChain::operator()(int parent, int child, const DependencyGraph& dg) const {
    if(! ft_.is_verb(child - 1)) return false;
    const std::string& edge = dg.edge(parent,child);
    if(edge == "aux" || edge =="auxpass") return false;
    return true;
  }
  const std::vector<std::string>& WordChain::get() const {
    return chain_;
  }
}
